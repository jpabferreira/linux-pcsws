#include "sched.h"
#include "cpudl.h"
#include <linux/slab.h>
#include <linux/syscalls.h>
#include <linux/math64.h>

inline void init_pcsws_job(struct rq *rq, struct task_struct *p) {
	RB_CLEAR_NODE(&p->pcsws.pjob_node);
    p->pcsws.pcsws_job.release = 0;
}

/* Returns 1 if @pcsws_se is the leftmost pjob in @pcsws_rq */
static inline int is_leftmost_pjob(struct pcsws_rq *pcsws_rq, struct sched_pcsws_entity *pcsws_se)
{
    return pcsws_rq->leftmost_pjob == &pcsws_se->pjob_node;
}

/* Returns 1 if @a < @b */
static inline int dl_time_before(u64 a, u64 b)
{
    return (s64)(a - b) < 0;
}

/*
 * Tells if entity @a should preempt entity @b. That's true when:
 *
 *  - Above conditions are met;
 *  - Entity @a is marked with SF_HEAD flag;
 *  - Entity @a has equal deadline than entity @b and @b is not a 'special task'.
 */
static inline int entity_preempt_equal_pcsws(struct sched_pcsws_entity *a, struct sched_pcsws_entity *b)
{
    return (a->pcsws_job.deadline == b->pcsws_job.deadline || dl_time_before(a->pcsws_job.deadline, b->pcsws_job.deadline));
}

/* Returns 1 if @pcsws_se is in the global rq */
static inline int on_global_rq(struct sched_pcsws_entity *pcsws_se)
{
    return !RB_EMPTY_NODE(&pcsws_se->rq_node);
}

/* Tasks have missed their previous deadline when their next release time is in the past */
inline int deadline_missed_pcsws(struct rq *rq, struct sched_pcsws_entity *pcsws_se)
{
    int dmiss;
    u64 damount;

    dmiss = dl_time_before(pcsws_se->pcsws_job.release, rq->clock);

    if (!dmiss)
        return 0;

    /*
     * Record statistics about maximum deadline
     * misses.
     */
    damount = rq->clock - pcsws_se->pcsws_job.release;

    //pcsws_se->stats.dmiss_max = max(pcsws_se->stats.dmiss_max, damount);
    pcsws_se->stats.dmiss_max = 1;

    return 1;
}

/* Returns the container struct task_struct of @pcsws_se */
static inline struct task_struct *task_of_pcsws_se(struct sched_pcsws_entity  *pcsws_se)
{
    if (pcsws_se)
        return container_of(pcsws_se, struct task_struct, pcsws);

    return NULL;
}

static inline struct rq *rq_of_pcsws_se(struct sched_pcsws_entity *pcsws_se)
{
    struct task_struct *p = task_of_pcsws_se(pcsws_se);
    struct rq *rq = task_rq(p);

    return rq;
}

/* Returns the container struct rq of @pcsws_rq */
static inline struct rq *rq_of_pcsws_rq(struct pcsws_rq *pcsws_rq)
{
    return container_of(pcsws_rq, struct rq, pcsws);
}

/* Returns 1 if @pcsws_se is in a local pcsws runqueue */
static inline int on_pcsws_rq(struct sched_pcsws_entity *pcsws_se)
{
    return !RB_EMPTY_NODE(&pcsws_se->pjob_node);
}

static int dispatch_pcsws(struct rq *rq, struct task_struct *p);

static void replenish_pcsws_entity(struct rq *rq, struct sched_pcsws_entity *pcsws_se)
{

    //BUG_ON(pi_se->pcsws_runtime <= 0);

    /*
     * This could be the case for a !-pcsws task that is boosted.
     * Just go with full inherited parameters.

    if (pcsws_se->pcsws_deadline == 0) {
        pcsws_se->deadline = rq->clock + pi_se->pcsws_deadline;
        pcsws_se->runtime = pi_se->pcsws_runtime;
    }*/

    /*
     * We keep moving the deadline away until we get some
     * available runtime for the entity. This ensures correct
     * handling of situations where the runtime overrun is
     * arbitrary large.
     */
    while (pcsws_se->pcsws_job.charge <= 0) {
        pcsws_se->pcsws_job.deadline += pcsws_se->pcsws_ded.period;
        pcsws_se->pcsws_job.charge += pcsws_se->pcsws_ded.budget;
    }

    /*
     * At this point, the deadline really should be "in
     * the future" with respect to rq->clock. If it's
     * not, we are, for some reason, lagging too much!
     * Anyway, after having warn userspace abut that,
     * we still try to keep the things running by
     * resetting the deadline and the budget of the
     * entity.
     */
    if (dl_time_before(pcsws_se->pcsws_job.deadline, rq->clock)) {
        pcsws_se->pcsws_job.deadline = rq->clock + pcsws_se->pcsws_ded.deadline;
        pcsws_se->pcsws_job.charge = pcsws_se->pcsws_ded.budget;
    }
}

/*
 * We are being explicitly informed that a new instance is starting,
 * and this means that:
 *  - the absolute deadline of the entity has to be placed at
 *    current time + relative deadline;
 *  - the runtime of the entity has to be set to the maximum value.
 *
 * The capability of specifying such event is useful whenever a -pcsws
 * entity wants to (try to!) synchronize its behaviour with the scheduler's
 * one, and to (try to!) reconcile itself with its own scheduling
 * parameters.
 */
static inline void setup_new_pcsws_entity(struct rq *rq, struct task_struct *p)
{
    WARN_ON(!p->pcsws.pcsws_new /*|| p->pcsws.throttled*/);

    //p->pcsws.stolen = -1;
    //p->pcsws.nr_pjobs = 0;
    p->pcsws.pcsws_new = 0;
    p->pcsws.stats.tot_runtime = 0;
    atomic_inc(&p->pcsws.pcsws_job.nr);

    replenish_pcsws_entity(rq, &p->pcsws);
}

static bool pcsws_entity_overflow(struct sched_pcsws_entity *pcsws_se, u64 t)
{
   u64 left, right;

   struct pcsws_job *job = &pcsws_se->pcsws_job;
   struct pcsws_dedicated *server = &pcsws_se->pcsws_ded;

   /*
    * left and right are the two sides of the equation above,
    * after a bit of shuffling to use multiplications instead
    * of divisions.
    *
    * Note that none of the time values involved in the two
    * multiplications are absolute: dl_deadline and dl_runtime
    * are the relative deadline and the maximum runtime of each
    * instance, runtime is the runtime left for the last instance
    * and (deadline - t), since t is rq->clock, is the time left
    * to the (absolute) deadline. Even if overflowing the u64 type
    * is very unlikely to occur in both cases, here we scale down
    * as we want to avoid that risk at all. Scaling down by 10
    * means that we reduce granularity to 1us. We are fine with it,
    * since this is only a true/false check and, anyway, thinking
    * of anything below microseconds resolution is actually fiction
    * (but still we want to give the user that illusion >;).
    */
   left = (server->period >> PCSWS_SCALE) * (job->charge >> PCSWS_SCALE);
   right = ((job->deadline - t) >> PCSWS_SCALE) *
       (server->budget >> PCSWS_SCALE);

   return dl_time_before(right, left);
}

/* Update @pcsws_se data to set up a new instance
 *
 * When a -pCSWS entity is queued back on the runqueue, its runtime and
 * deadline might need updating.
 *
 * The policy here is that we update the deadline of the entity only if:
 *  - the current deadline is in the past,
 *  - using the remaining runtime with the current deadline would make
 *    the entity exceed its bandwidth.
 */
inline void update_task_pcsws(struct rq *rq, struct sched_pcsws_entity *pcsws_se)
{
    struct task_struct *p = task_of_pcsws_se(pcsws_se);

    /*
     * The arrival of a new instance needs special treatment, i.e.,
     * the actual scheduling parameters have to be "renewed".
     */
    if (pcsws_se->pcsws_new) {
      //-  printk(KERN_EMERG "Setup new\n");
        setup_new_pcsws_entity(rq, p);
        return;
    }

    if (dl_time_before(pcsws_se->pcsws_job.deadline, rq->clock) ||
        pcsws_entity_overflow(pcsws_se, rq->clock)) {
        pcsws_se->pcsws_job.deadline = rq->clock + pcsws_se->pcsws_ded.deadline;
        pcsws_se->pcsws_job.charge = pcsws_se->pcsws_ded.budget;
    }
}

/* Returns the leftmost pjob in @pcsws_rq */
static struct sched_pcsws_entity *pick_next_pjob_pcsws(struct pcsws_rq *pcsws_rq)
{
    struct rb_node *edf = pcsws_rq->leftmost_pjob;

    if (!edf)
        return NULL;

    return rb_entry(edf, struct sched_pcsws_entity, pjob_node);
}

/* Returns the highest priority task in the global runqueue */
static struct sched_pcsws_entity *__pick_next_task_pcsws(struct pcsws_ready *global_rq)
{
    struct rb_node *edf = global_rq->leftmost;

    if (!edf)
        return NULL;

    return rb_entry(edf, struct sched_pcsws_entity, rq_node);
}

/* Enqueue @p into the global runqueue */
static void __enqueue_task_pcsws(struct pcsws_ready *global_rq, struct sched_pcsws_entity *pcsws_se)
{
    struct rb_node **link = &global_rq->tasks.rb_node;
    struct rb_node *parent = NULL;
    struct sched_pcsws_entity *entry;
    int edf = 1;

    //BUG_ON(on_global_rq(pcsws_se));

    /* TODO: Optimization in case enqueueing task is next edf */
    while (*link) {
        parent = *link;
        entry = rb_entry(parent, struct sched_pcsws_entity, rq_node);

        if (dl_time_before(pcsws_se->pcsws_job.deadline, entry->pcsws_job.deadline))
            link = &parent->rb_left;
        else {
            link = &parent->rb_right;
            edf = 0;
        }
    }

    if (edf)
        global_rq->leftmost = &pcsws_se->rq_node;

    rb_link_node(&pcsws_se->rq_node, parent, link);
    rb_insert_color(&pcsws_se->rq_node, &global_rq->tasks);

    global_rq->nr_running++;
}

/* Push @p into the global runqueue */
static void push_task_pcsws(struct rq *rq, struct task_struct *p, int preempted)
{
    struct pcsws_ready *global_rq = rq->pcsws.pcsws_ready;

    if (preempted){
      //-  printk(KERN_EMERG "[Preempted]\n");
        deactivate_task(rq, p, 0);
    }

    raw_spin_lock(&global_rq->lock);
    __enqueue_task_pcsws(global_rq, &p->pcsws);
    raw_spin_unlock(&global_rq->lock);
}

/*
 * Tries to push a -pcsws task to a "random" idle rq.
 */
static int push_idle_pcsws(struct rq *this_rq, struct task_struct *p)
{
    struct rq *target_rq;
    int ret = 0, target_cpu;
    struct cpudl *cp = &this_rq->rd->pcswsc_cpudl;

retry:
    target_cpu = find_idle_cpu_pcsws(cp);

    if (target_cpu == -1)
        return 0;

    //printk(KERN_INFO "idle cpu %d\n", target_cpu);

    target_rq = cpu_rq(target_cpu);

    /* We might release rq lock */
    get_task_struct(p);

    double_lock_balance(this_rq, target_rq);

    if (unlikely(target_rq->pcsws.nr_running)) {
        double_unlock_balance(this_rq, target_rq);
        put_task_struct(p);
        target_rq = NULL;
        goto retry;
    }

    /*if (p->pcsws.parent || p->pcsws.stolen != -1) {
        //printk(KERN_INFO "--steal\n");
        //deactivate_task(this_rq, p, 0);
        if (p->state != TASK_WAKING)
            p->pcsws.stolen = this_rq->cpu;
        target_rq->pcsws.tot_steals++;
    }*/

    set_task_cpu(p, target_cpu);
    activate_task(target_rq, p, 0);

    ret = 1;
    resched_task(target_rq->curr);

    double_unlock_balance(this_rq, target_rq);
    put_task_struct(p);

    return ret;
}

/* Push @p into the lowest priority CPU's runqueue */
static int push_latest_pcsws(struct rq *this_rq, struct task_struct *p, int target_cpu)
{
    struct rq *target_rq;
    int ret = 0;

    target_rq = cpu_rq(target_cpu);

    /* We might release rq lock */
    get_task_struct(p);

    double_lock_balance(this_rq, target_rq);

    /* TODO check if in the meanwhile a task was dispatched to here */
    if (target_rq->pcsws.nr_running && !dl_time_before(p->pcsws.pcsws_job.deadline, target_rq->pcsws.earliest_dl)) {
        double_unlock_balance(this_rq, target_rq);
        put_task_struct(p);

        return ret;
    }

    deactivate_task(this_rq, p, 0);
    set_task_cpu(p, target_cpu);
    activate_task(target_rq, p, ENQUEUE_HEAD);
    resched_task(target_rq->curr);

    double_unlock_balance(this_rq, target_rq);




    put_task_struct(p);


    ret = 1;



    return ret;
}

/* Place the newly awakened task @p in the appropriate runqueue */
static int dispatch_pcsws(struct rq *rq, struct task_struct *p)
{
    struct sched_pcsws_entity *edf;
    int target_cpu;
    struct cpudl *cp = &rq->rd->pcswsc_cpudl;





    goto global;






    edf = __pick_next_task_pcsws(rq->pcsws.pcsws_ready);

    if (edf && !dl_time_before(p->pcsws.pcsws_job.deadline, edf->pcsws_job.deadline))
        goto global;
    /* else -> @p is EDF */

    target_cpu = find_latest_cpu_pcsws(cp);

    /* If @p is to stay allocated to the current CPU */
    if (target_cpu == rq->cpu) {
        /* If @p is not the highest priority task in the runqueue */
        if (rq->pcsws.nr_running && !dl_time_before(p->pcsws.pcsws_job.deadline, rq->pcsws.earliest_dl))
            goto try_latest;

        //printk(KERN_EMERG "Don't push\n");
        return 0;
    }

try_latest:
    if (target_cpu == -1)
        goto global;

    /* Task p is to be pushed to the least priority CPU */
    if (push_latest_pcsws(rq, p, target_cpu)) {
        //printk(KERN_EMERG "Push Latest\n");
        return 1;
    }

global:
    /* Task p is to be pushed to the global runqueue */
    push_task_pcsws(rq, p, 0);
    //printk(KERN_EMERG "Push global\n");

    return 1;
}

/* Start tasks' timer, setting the next release time */
static inline int start_timer_pcsws(struct sched_pcsws_entity *pcsws_se)
{
    struct rq *rq = rq_of_pcsws_se(pcsws_se);

    ktime_t now, act;
    ktime_t soft, hard;
    unsigned long range;
    s64 delta;

    /*if (boosted)
        return 0;
    */
    /*
     * We want the timer to fire at the deadline, but considering
     * that it is actually coming from rq->clock and not from
     * hrtimer's time base reading.
     */
    act = ns_to_ktime(pcsws_se->pcsws_job.deadline);
    now = hrtimer_cb_get_time(&pcsws_se->timer);
    delta = ktime_to_ns(now) - rq->clock;
    act = ktime_add_ns(act, delta);

    //printk(KERN_EMERG "[Tryig o set timer]");

    /*
     * If the expiry time already passed, e.g., because the value
     * chosen as the deadline is too small, don't even try to
     * start the timer in the past!
     */
    if (ktime_us_delta(act, now) < 0)
        return 0;

    hrtimer_set_expires(&pcsws_se->timer, act);

    soft = hrtimer_get_softexpires(&pcsws_se->timer);
    hard = hrtimer_get_expires(&pcsws_se->timer);
    range = ktime_to_ns(ktime_sub(hard, soft));
    __hrtimer_start_range_ns(&pcsws_se->timer, soft,
                 range, HRTIMER_MODE_ABS, 0);

  //-  printk(KERN_EMERG "[Timer set to] %llu\n", ktime_to_ns(act));

    return hrtimer_active(&pcsws_se->timer);
}

/* Update @pcsws_rq and CPU heap data when @pcsws_se is dequeued */
static
void dec_pjobs_pcsws(struct pcsws_rq *pcsws_rq, struct sched_pcsws_entity *pcsws_se)
{
    struct rq *rq = rq_of_pcsws_rq(pcsws_rq);
    struct task_struct *p = task_of_pcsws_se(pcsws_se);
    struct sched_pcsws_entity *stealable;
    struct pcsws_ready * global_rq = pcsws_rq->pcsws_ready;
    int ret = 0;

    WARN_ON(!pcsws_prio(p->prio));
    WARN_ON(!pcsws_rq->nr_running);
    pcsws_rq->nr_running--;

    if (!pcsws_rq->nr_running) {
        /* If there are no more pjobs to run, we declare this rq idle  */
        pcsws_rq->earliest_dl = 0;
        cpudl_set(&rq->rd->pcswsc_cpudl, rq->cpu, 0, 0);
        smp_wmb();
    } else {
        /*
        if (!RB_EMPTY_NODE(&pcsws_se->stealable_pjob_node))
            return;
        */
        if (pcsws_rq->earliest_dl != pcsws_se->pcsws_job.deadline)
            return;
        /*
        if (!has_stealable_pjobs(pcsws_rq))
            return;
        */

        /* The leftmost stealable pjob and our next pjob share the same deadline
        stealable = rb_entry(pcsws_rq->leftmost_stealable_pjob, struct sched_pcsws_entity, stealable_pjob_node);

        if (stealable->job.deadline > pcsws_rq->earliest_dl) {
            /*
             * If the next pjob has lower priority than the highest priority task on
             * global rq, we try to pull that task.
             * Clearing the earlieast deadline value (earliest_dl=0) on the rq is
             * a trick to allow next's rq statistics update, keeping the current ones.

            if (priority_inversion_pcsws(pcsws_rq, stealable)) {
                pcsws_rq->earliest_dl = 0;

                raw_spin_lock(&global_rq->lock);
                ret = pull_task_pcsws(rq);
                raw_spin_unlock(&global_rq->lock);

                if (ret)
                    return;
            }
        }
        */

        /*
         * We update statistics about this rq when next
         * pjob has a different deadline than the dequeueing one.

        if (pcsws_rq->earliest_dl != stealable->job.deadline) {
            pcsws_rq->earliest_dl = stealable->job.deadline;
            cpupcsws_set(&rq->rd->pcswsc_cpudl, rq->cpu, pcsws_rq->earliest_dl, 0, 1);
            smp_wmb();
        }*/
    }
}

/* Dequeue @pcsws_se from the global runqueue @global_rq */
static void __dequeue_task_pcsws(struct pcsws_ready *global_rq, struct sched_pcsws_entity *pcsws_se)
{
    if (!on_global_rq(pcsws_se))
        return;

    global_rq->nr_running--;

    if (global_rq->leftmost == &pcsws_se->rq_node) {
        struct rb_node *next_node;

        next_node = rb_next(&pcsws_se->rq_node);
        global_rq->leftmost = next_node;
    }

    rb_erase(&pcsws_se->rq_node, &global_rq->tasks);
    RB_CLEAR_NODE(&pcsws_se->rq_node);
}

/* Dequeue @pcsws_se from the local runqueue @pcsws_rq */
static void __dequeue_pjob_pcsws(struct pcsws_rq *pcsws_rq, struct sched_pcsws_entity *pcsws_se)
{
    if (!on_pcsws_rq(pcsws_se))
        return;

    if (pcsws_rq->leftmost_pjob == &pcsws_se->pjob_node) {
        struct rb_node *next_node;

        next_node = rb_next(&pcsws_se->pjob_node);
        pcsws_rq->leftmost_pjob = next_node;
    }

    rb_erase(&pcsws_se->pjob_node, &pcsws_rq->pjobs);
    RB_CLEAR_NODE(&pcsws_se->pjob_node);

    dec_pjobs_pcsws(pcsws_rq, pcsws_se);
}

/* Dequeue pjob @p */
static void dequeue_pjob_pcsws(struct rq *rq, struct task_struct *p)
{
    struct sched_pcsws_entity *pcsws_se = &p->pcsws;

    __dequeue_pjob_pcsws(&rq->pcsws, pcsws_se);
}

/* Scheduling class function  (dequeue_task_pcsws) */
static void dequeue_task_pcsws(struct rq *rq, struct task_struct *p, int flags)
{
    dequeue_pjob_pcsws(rq, p);
}

/* Move the highest priority task from the global runqueue to the local runqueue @this_rq */
static int pull_task_pcsws(struct rq *this_rq)
{
    struct task_struct *p;
    struct sched_pcsws_entity *pcsws_se;
    struct pcsws_ready *global_rq = this_rq->pcsws.pcsws_ready;


    if (!global_rq->nr_running)
        return 0;

    pcsws_se = __pick_next_task_pcsws(global_rq);

    if (unlikely(!pcsws_se))
        return 0;

    p = task_of_pcsws_se(pcsws_se);
    WARN_ON(!pcsws_task(p));

    /* We need to transfer the task from the global runqueue to this runqueue */
    __dequeue_task_pcsws(global_rq, pcsws_se);

    /*
    if (pcsws_se->stolen == this_rq->cpu)
        pcsws_se->stolen = -1;*/

    set_task_cpu(p, this_rq->cpu);
    activate_task(this_rq, p, ENQUEUE_HEAD);

    //this_rq->pcsws.tot_pulls++;

    return 1;
}



/* Scheduling class function (select_task_rq) */
static int select_task_rq_pcsws(struct task_struct *p, int sd_flag, int flags)
{
    unsigned long rq_flags;
    struct rq *rq = task_rq_lock(p, &rq_flags);
    struct task_struct *curr = rq->curr;
    struct cpudl *cp = &rq->rd->pcswsc_cpudl;
    int cpu, target_cpu;

    cpu = task_cpu(p);

    if (unlikely(!pcsws_task(p)))
        goto out;

    if (on_pcsws_rq(&p->pcsws))
        goto out;

    if (sd_flag != SD_BALANCE_WAKE && sd_flag != SD_BALANCE_FORK)
        goto out;

    /* If we're not dealing with a -pcsws task, we don't even bother */
    if (unlikely(!pcsws_task(p))) {
        goto out;
    }

    /*
     * If the current CPU is idle at this point, there's no need to migrate.
     * Notice that we consider the CPU idle whenever it is not executing -pcsws
     * tasks. As far as our scheduling class is concerned, whether a CPU is
     * idle or executing tasks from other scheduling classes, makes no real
     * difference
     */
    if (!rq->pcsws.nr_running)
        goto out;

    /*
     * Having failed to wake up the task in this rq, we'll try to find an idle
     * processor for it.
     */
    target_cpu = find_idle_cpu_pcsws(cp);
    if (target_cpu != -1) {
        cpu = target_cpu;
        goto out;
    }

    /*
     * We can't find an idle cpu, so we'll try to place the task on the lowest
     * priority rq.
     */
    target_cpu = find_latest_cpu_pcsws(cp);
    if (target_cpu != -1)
        cpu = target_cpu;

    /*
     * If all CPUs are busy and our task will not cause preemption, we'll wake it
     * up in the current CPU and push it to the global rq when enqueue_task_pcsws
     * is called.
     */
out:
    task_rq_unlock(rq, p, &rq_flags);

	printk(KERN_EMERG "[select_task_rq] %d\n", cpu);

    return cpu;
}

/*
 * Scheduling class function (check_preempt_curr)
 * Check if @p is to be preempted
 */
static void check_preempt_curr_pcsws(struct rq *rq, struct task_struct *p, int flags)
{
    if (!is_leftmost_pjob(&rq->pcsws, &rq->curr->pcsws))
            resched_task(rq->curr);
}

/* Update @pcsws_rq and CPU heap data when @pcsws_se is enqueued */
static
void inc_pjobs_pcsws(struct pcsws_rq *pcsws_rq, struct sched_pcsws_entity *pcsws_se)
{
    struct rq *rq = rq_of_pcsws_rq(pcsws_rq);
    struct task_struct *p = task_of_pcsws_se(pcsws_se);
    u64 deadline = pcsws_se->pcsws_job.deadline;

    WARN_ON(!pcsws_prio(p->prio));
    pcsws_rq->nr_running++;

    /*
     * We update statistics about this rq only when
     * enqueueing pjob has lower deadline than our current
     * one, if we got any running on, that's it.
     */
    if (pcsws_rq->earliest_dl == 0 || dl_time_before(deadline, pcsws_rq->earliest_dl)) {
        pcsws_rq->earliest_dl = deadline;

        cpudl_set(&rq->rd->pcswsc_cpudl, rq->cpu, deadline, 1);
        smp_wmb();
    }
}

/*
 * Enqueue the pjob @pcsws_se into the local runqueue @pcsws_rq
 * The list of ready -pcsws pjobs is a rb-tree ordered by deadline.
 * Pjobs with equal deadline obey to a LIFO order.
 */
static void __enqueue_pjob_pcsws(struct pcsws_rq *pcsws_rq, struct sched_pcsws_entity *pcsws_se)
{
    struct rb_node **link = &pcsws_rq->pjobs.rb_node;
    struct rb_node *parent = NULL;
    struct sched_pcsws_entity *entry;
    int edf = 1;

  //-  printk(KERN_EMERG "[Trying to enqueue]\n");
    //BUG_ON(on_pcsws_rq(pcsws_se));
    if (on_pcsws_rq((pcsws_se))) {
        printk(KERN_EMERG "[Bug: On rq]\n");
        return;
    }

    /* TODO: Optimization in case enqueueing task is next edf */


    while (*link) {
        parent = *link;
        entry = rb_entry(parent, struct sched_pcsws_entity, pjob_node);

        if (/*(pcsws_se->help_first && !is_leftmost_pjob(pcsws_rq, entry) && entity_preempt_equal_pcsws(pcsws_se, entry)) ||
            (pcsws_se->help_first && is_leftmost_pjob(pcsws_rq, entry) && entity_preempt_pcsws(pcsws_se, entry)) ||
            (!pcsws_se->help_first && */
            entity_preempt_equal_pcsws(pcsws_se, entry))/*)*/
            link = &parent->rb_left;
        else {
            link = &parent->rb_right;
            edf = 0;
        }
    }

    if (edf)
        pcsws_rq->leftmost_pjob = &pcsws_se->pjob_node;

    rb_link_node(&pcsws_se->pjob_node, parent, link);
    rb_insert_color(&pcsws_se->pjob_node, &pcsws_rq->pjobs);

    //pcsws_rq->tot_enqueues++;
    inc_pjobs_pcsws(pcsws_rq, pcsws_se);

  //-  printk(KERN_EMERG "[Enqueued]\n");
}

/* Enqueue pjob @p */
static void enqueue_pjob_pcsws(struct rq *rq, struct task_struct *p, int flags)
{
    struct sched_pcsws_entity *pcsws_se = &p->pcsws;
    int new_job = pcsws_se->pcsws_new;

    /*
     * If this is a wakeup or a new instance, the scheduling
     * parameters of the task might need updating. Otherwise,
     * we want a replenishment of its runtime.
     */
    if (!new_job && flags & ENQUEUE_REPLENISH)
        replenish_pcsws_entity(rq, pcsws_se);
    else
        update_task_pcsws(rq, pcsws_se);

    __enqueue_pjob_pcsws(&rq->pcsws, pcsws_se);
    /*
    if (!task_current(rq, p) && rq->pcsws.nr_running > 1 && !(flags & ENQUEUE_HEAD))
        enqueue_stealable_pjob_pcsws(&rq->pcsws, pcsws_se); */
}

/* Scheduling class function (enqueue_task_pcsws) */
static void enqueue_task_pcsws(struct rq *rq, struct task_struct *p, int flags)
{
    struct pcsws_rq *pcsws_rq = &rq->pcsws;
    struct pcsws_ready *global_rq = pcsws_rq->pcsws_ready;
    //printk(KERN_EMERG "[Throttled] %d\n", p->pcsws.throttled);
    /*
     * If p is throttled, we do nothing. In fact, if it exhausted
     * its budget it needs a replenishment and, since it now is on
     * its rq, the bandwidth timer callback (which clearly has not
     * run yet) will take care of this.
     */
    if (p->pcsws.throttled)
        return;

    /*
     * Provided the current CPU is idle, we know the task will be enqueued
     */
    if (!rq->pcsws.nr_running)
        goto enqueue_local;

    /*
     * The current CPU is not idle, so we'll only enqueue the task if
     * rq->curr is to be preempted
     */
    if (!entity_preempt_equal_pcsws(&p->pcsws, &rq->curr->pcsws)) {
		printk(KERN_EMERG "[Enqueue global]\n");
        __enqueue_task_pcsws(global_rq, &p->pcsws);
        return;
    }


enqueue_local:
	printk(KERN_EMERG "[Enqueue local]\n");
    enqueue_pjob_pcsws(rq, p, flags);

    check_preempt_curr_pcsws(rq, rq->curr, 0);
}

/* Tasks' timer function */
static enum hrtimer_restart timer_pcsws(struct hrtimer *timer)
{
    unsigned long flags;
    struct sched_pcsws_entity *pcsws_se = container_of(timer,
                             struct sched_pcsws_entity,
                             timer);

    struct task_struct *p = task_of_pcsws_se(pcsws_se);
    struct rq *rq = task_rq_lock(p, &flags);

  //-  printk(KERN_EMERG "Timer fired\n");

    /*
     * We need to take care of a possible races here. In fact, the
     * task might have changed its scheduling policy to something
     * different from SCHED_PCSWS (through sched_setscheduler()).
     */
    /*
    if (!pcsws_task(p) || pcsws_se->pcsws_new) {
      //-  printk(KERN_EMERG "Timer did nothing 1\n");
        goto unlock;
    }*/
    if (!pcsws_task(p)) {
      //-  printk(KERN_EMERG "Timer did nothing 1\n");
        goto unlock;
    }
    if (pcsws_se->pcsws_new) {
      //-  printk(KERN_EMERG "Timer did nothing 2\n");
        goto unlock;
    }

    pcsws_se->throttled = 0;
    if (p->on_rq) {
        //printk(KERN_EMERG "On RQ\n");
        enqueue_task_pcsws(rq, p, ENQUEUE_REPLENISH);
        if (task_has_pcsws_policy(rq->curr)) {
          //-  printk(KERN_EMERG "Before check_preemt\n");
            check_preempt_curr_pcsws(rq, p, 0);
          //-  printk(KERN_EMERG "After check_preemt");
        }
        else {
            resched_task(rq->curr);
        }
    }
unlock:
    task_rq_unlock(rq, p, &flags);

    return HRTIMER_NORESTART;
}

/* Init tasks' timer */
void init_timer_pcsws(struct sched_pcsws_entity *pcsws_se)
{
    struct hrtimer *timer = &pcsws_se->timer;

    if (hrtimer_active(timer)) {
        hrtimer_try_to_cancel(timer);
        return;
    }

    hrtimer_init(timer, CLOCK_MONOTONIC, HRTIMER_MODE_REL);
    timer->function = timer_pcsws;
}

/* Scheduling class function (pick_next_task_pcsws)
 * Returns the highest priority task in the local runqueue
 * If the local runqueue contains no running pjobs, we pull the highest priority pjob
 * from the global runqueue
 */
static struct task_struct *pick_next_task_pcsws(struct rq *rq)
{
    struct sched_pcsws_entity *pcsws_se;
    struct task_struct *next;
    struct pcsws_rq *pcsws_rq = &rq->pcsws;
    struct pcsws_ready *global_rq = pcsws_rq->pcsws_ready;
    int ret;

    /* TODO: if new pjobs spawn after we get a null here, we wont be aware of them until next pick round */
    if (!pcsws_rq->nr_running) {
        raw_spin_lock(&global_rq->lock);
        ret = pull_task_pcsws(rq);
        raw_spin_unlock(&global_rq->lock);

        if (!ret)
            /*if (!steal_pjob_pcsws(rq))*/
                return NULL;
        pcsws_se = pick_next_pjob_pcsws(pcsws_rq);
        //printk(KERN_EMERG "[Pulled] %d\n", atomic_read(&pcsws_se->pcsws_job.nr));
    }

    pcsws_se = pick_next_pjob_pcsws(pcsws_rq);
    BUG_ON(!pcsws_se);

    next = task_of_pcsws_se(pcsws_se);
    BUG_ON(!pcsws_task(next));

    next->se.exec_start = rq->clock_task;

    if (next) {
      //-  printk(KERN_EMERG "[Returned] %d\n", atomic_read(&pcsws_se->pcsws_job.nr));
    }

    return next;
}

static
int dl_runtime_exceeded(struct rq *rq, struct sched_pcsws_entity *pcsws_se)
{
    int dmiss = dl_time_before(pcsws_se->pcsws_job.deadline, rq->clock);
    int rorun = pcsws_se->pcsws_job.charge <= 0;

    if (/*pcsws_se->pcsws_new == 1 ||*/ (!rorun && !dmiss))
        return 0;



    /*
     * Record statistics about last and maximum deadline
     * misses and runtime overruns.
     */
    /*
    if (dmiss) {
        u64 damount = rq->clock - pcsws_se->deadline;

        schedstat_set(pcsws_se->stats.last_dmiss, damount);
        schedstat_set(pcsws_se->stats.dmiss_max,
                  max(pcsws_se->stats.dmiss_max, damount));
    }
    if (rorun) {
        u64 ramount = -pcsws_se->runtime;

        schedstat_set(pcsws_se->stats.last_rorun, ramount);
        schedstat_set(pcsws_se->stats.rorun_max,
                  max(pcsws_se->stats.rorun_max, ramount));
    }
    */

    /*
     * If we are beyond our current deadline and we are still
     * executing, then we have already used some of the runtime of
     * the next instance. Thus, if we do not account that, we are
     * stealing bandwidth from the system at each deadline miss!
     */
    if (dmiss) {
        pcsws_se->pcsws_job.charge = rorun ? pcsws_se->pcsws_job.charge : 0;
        pcsws_se->pcsws_job.charge -= rq->clock - pcsws_se->pcsws_job.deadline;
    }

    return 1;
}

/* Update current task data and statistics */
static void update_curr_pcsws(struct rq *rq)
{
    struct task_struct *curr = rq->curr;
    struct sched_pcsws_entity *pcsws_se = &curr->pcsws;
    u64 delta_exec;


    //printk(KERN_EMERG "[+");
    if (!pcsws_task(curr) || !on_pcsws_rq(pcsws_se))
        return;
    //printk(KERN_EMERG "-]");

    delta_exec = rq->clock_task - curr->se.exec_start;
    if (unlikely((s64)delta_exec < 0))
        delta_exec = 0;

    /*schedstat_set(curr->se.statistics.exec_max,
              max(curr->se.statistics.exec_max, delta_exec));
    */

    curr->se.sum_exec_runtime += delta_exec;
    //schedstat_add(&rq->pcsws, exec_clock, delta_exec);
    /* Maintain exec runtime for a thread group. @sched_stats.h */
    account_group_exec_runtime(curr, delta_exec);

    curr->se.exec_start = rq->clock_task;
    /* Charge this task's execution time to its accounting group. @sched.c */
    cpuacct_charge(curr, delta_exec);
    sched_rt_avg_update(rq, delta_exec);

    //pcsws_se->stats.tot_runtime += delta_exec;

    pcsws_se->pcsws_job.charge -= delta_exec;


    /*
     * In case we arrive here following a call to wait_interval_pcsws, having
     * accounted for the time elapsed since the last tick, we don't need to
     * test if the task has exceeded its timing constraints.
     * The task will yield the processor and wake up at the next release time,
     * where new scheduling parameters will be assigned.
     */
    if (pcsws_se->pcsws_new)
        return;

    /*
     * If pcsws_se has exceeded its runtime we take it out of the runqueue,
     * activate the throttled flag and try to set the timer which, upon firing,
     * will update the scheduling parameters and enqueue it back.
     * In case we don't succeed at setting the timer we know that the thread
     * has missed its deadline, so it's necessary to update its scheduling
     * parameters by replenishing the budget and postponing the deadline.
     */
    if (dl_runtime_exceeded(rq, pcsws_se)) {
      //-  printk(KERN_EMERG "[dl_runtime_exceeded]    [Clock] %llu     [Dl] %llu     [Charge]     %lld\n", rq->clock, curr->pcsws.pcsws_job.deadline, curr->pcsws.pcsws_job.charge);
        __dequeue_pjob_pcsws(&rq->pcsws, &curr->pcsws);
        if (start_timer_pcsws(pcsws_se/*, curr->pcsws.pcsws_boosted*/)) {
            pcsws_se->throttled = 1;
        }
        else {
            //printk(KERN_EMERG "[Enqueue replenish] %d", curr->pid);
            enqueue_task_pcsws(rq, curr, ENQUEUE_REPLENISH);
        }

        if (!is_leftmost_pjob(&rq->pcsws, pcsws_se)) {
          //-  printk(KERN_EMERG "[Did resched] %d\n", curr->pid);
            resched_task(curr);
        }
    }
}

/* Scheduling class function (put_prev_task_pcsws) */
static void put_prev_task_pcsws(struct rq *rq, struct task_struct *p)
{
    struct sched_pcsws_entity *pcsws_se = &p->pcsws;
    struct pcsws_ready *global_rq = rq->pcsws.pcsws_ready;
    int migrate = 0;

    update_curr_pcsws(rq);
    //p->se.exec_start = 0;

    /* A stolen pjob can not be delayed by a preemption like any other local pjob,
     * because that blocking time could make it lose its deadline. That is, instead of helping
     * it to finish earlier, this rq would add a substancial response time, which could lead
     * in a execution time higher than WCET (sequential-based): a deadline violation certainly
     * avoidable, if the pjob hasn't been stolen.
     * By forcing pjob=0, we make sure this pjob ends up on global rq, thus becoming able
     * to preempt other cores or to be selected on the 2nd highest priority picking up phase.
     * */

    /*if (pcsws_se->stolen != -1)
        migrate = 1;
    */

    if (on_pcsws_rq(pcsws_se) && rq->pcsws.nr_running > 1) {
        if ((/*!pcsws_se->parent &&*/ pcsws_se->nr_pjobs == 0) /*|| migrate*/) {
          //-  printk(KERN_EMERG "Push task no put_prev");
            push_task_pcsws(rq, p, 1);
        }
        /*else
            enqueue_stealable_pjob_pcsws(&rq->pcsws, pcsws_se); */
    }

    /*if (!pcsws_se->parent && p->state == TASK_DEAD) {
        //printk(KERN_INFO "Hrtimer cancel 1\n");
        hrtimer_try_to_cancel(&pcsws_se->timer);
        if (pcsws_se->stats.dmiss_max > 0) {
            if (global_rq->random)
                trace_sched_task_stat_pcsws_random(prev, rq->clock);
            else
                trace_sched_task_stat_pcsws_earliest(prev, rq->clock);
        }
    }*/
}

/* Scheduling class function (yield_task_pcsws) */
static void yield_task_pcsws(struct rq *rq)
{
    struct task_struct *p = rq->curr;

    /*
     * We make the task go to sleep until its current deadline by
     * forcing its runtime to zero. This way, update_curr_pcsws() stops
     * it and the bandwidth timer will wake it up and will give it
     * new scheduling parameters (thanks to pcsws_new=1).
     */
    if (p->pcsws.pcsws_job.charge > 0) {
        rq->curr->pcsws.pcsws_new = 1;
        p->pcsws.pcsws_job.charge = 0;
    }

    update_curr_pcsws(rq);
}

/* Scheduling class function (set_cpud_allowed_pcsws) */
static void set_cpus_allowed_pcsws(struct task_struct *p,
                const struct cpumask *new_mask)
{
}

/* Scheduling class function (rq_online_pcsws) */
static void rq_online_pcsws(struct rq *rq)
{
}

/* Scheduling class function (rq_offline_pcsws) */
static void rq_offline_pcsws(struct rq *rq)
{
}

/* Scheduling class function (pre_schedule_pcsws) */
static void pre_schedule_pcsws(struct rq *rq, struct task_struct *prev)
{
}

/* Scheduling class function (post_schedule_pcsws) */
static void post_schedule_pcsws(struct rq *rq)
{
}

/* Scheduling class function (task_woken_pcsws) */
static void task_woken_pcsws(struct rq *rq, struct task_struct *p)
{
}

/* Scheduling class function (set_curr_task_pcsws) */
static void set_curr_task_pcsws(struct rq *rq)
{
    struct task_struct *p = rq->curr;

    update_task_pcsws(rq, &p->pcsws);
    p->se.exec_start = rq->clock_task;
}

/* Scheduling class function (task_tick_pcsws) */
static void task_tick_pcsws(struct rq *rq, struct task_struct *p, int queued)
{
    //update_curr_pcsws(rq);

    /*
#ifdef CONFIG_SCHED_HRTICK
    if (hrtick_enabled(rq) && queued && p->pcsws.pcsws_job.runtime > 0)
        start_hrtick_pcsws(rq, p);
#endif
*/
}

/* Scheduling class function (get_rr_interval_pcsws) */
static unsigned int get_rr_interval_pcsws(struct rq *rq, struct task_struct *task)
{
    return 0;
}

/* Scheduling class function (prio_changed_pcsws) */
static void
prio_changed_pcsws(struct rq *rq, struct task_struct *p, int oldprio)
{
}

/* Scheduling class function (switched_from_pcsws) */
static void switched_from_pcsws(struct rq *rq, struct task_struct *p)
{
    if (hrtimer_active(&p->pcsws.timer) && !pcsws_policy(p->policy)) {
        hrtimer_try_to_cancel(&p->pcsws.timer);
    }
}

/* Scheduling class function (switched_to_pcsws) */
static void switched_to_pcsws(struct rq *rq, struct task_struct *p)
{
    /*
     * If p is throttled, don't consider the possibility
     * of preempting rq->curr, the check will be done right
     * after its runtime will get replenished.
     */
    if (unlikely(p->pcsws.throttled))
        return;

    if (!rq->nr_running)
        check_preempt_curr_pcsws(rq, p, 0);
}

/* Scheduling class function (task_fork_pcsws) */
static void task_fork_pcsws(struct task_struct *p)
{
}

/* Scheduling class function (task_dead_pcsws) */
static void task_dead_pcsws(struct task_struct *p)
{
    struct sched_pcsws_entity *pcsws_se = &p->pcsws;
    hrtimer_try_to_cancel(&pcsws_se->timer);
}

/*
 * Here we check if --at time t-- a task (which is probably being
 * [re]activated or, in general, enqueued) can use its remaining runtime
 * and its current deadline _without_ exceeding the bandwidth it is
 * assigned (function returns true if it can).
 *
 * For this to hold, we must check if:
 *   pcsws_job.charge / (pcsws_job.deadline - t) < pcsws.budget / pcsws.deadline .
 */
static bool pcsws_check_bandwidth(struct sched_pcsws_entity *pcsws_se, u64 t)
{
    u64 left, right;

    /*
     * left and right are the two sides of the equation above,
     * after a bit of shuffling to use multiplications instead
     * of divisions.
     */
    left = pcsws_se->pcsws_ded.deadline * pcsws_se->pcsws_job.charge;
    right = (pcsws_se->pcsws_job.deadline - t) * pcsws_se->pcsws_ded.budget;

    return dl_time_before(left, right);
}

/*
 * This function makes the task sleep until at least the absolute time
 * instant specified in @rqtp.
 * The _at_least_ part comes from the fact that we want to be able
 * to give the task --as soon as it wakes-up-- its full runtime.
 * Therefore, if e.g. it is in overrun when this function is invoked,
 * of if @rqtp is too early, the sleeping time might be longer than asked.
 * It is intended to be used at the end of a periodic -deadline task
 * instance.
 */
long wait_interval_pcsws(struct task_struct *p, struct timespec *rqtp,
              struct timespec *rmtp)
{
    unsigned long flags;
    struct sched_pcsws_entity *pcsws_se = &p->pcsws;
    struct rq *rq = task_rq_lock(p, &flags);
    struct timespec lrqtp;
    u64 wakeup;

  //-  printk(KERN_EMERG "[Wait interval]     [Now] %llu     [Dl] %llu     [Charge]     %lld\n", rq->clock, p->pcsws.pcsws_job.deadline, p->pcsws.pcsws_job.charge);

    p->pcsws.pcsws_new = 1;
    update_curr_pcsws(rq);

    if (hrtimer_active(&pcsws_se->timer)) {
      //-  printk(KERN_EMERG "Timer was active\n");
        hrtimer_try_to_cancel(&pcsws_se->timer);
    }

    __dequeue_pjob_pcsws(&rq->pcsws, pcsws_se);

    /*
     * Task is asking for a new instance with full runtime but, since
     * it is in overrun, the only thing we can do is putting it to sleep
     * until the time it will have payed back for that (which could
     * be its next deadline or farther).

    if (pcsws_se->pcsws_job.charge < 0) {
        u64 runtime_exec = pcsws_se->pcsws.budget - pcsws_se->pcsws_job.charge;
        u64 rorun_ratio = runtime_exec / pcsws_se->pcsws.budget;

        WARN_ON(rorun_ratio == 0);

        wakeup = pcsws_se->pcsws_job.deadline + pcsws_se->pcsws.deadline * rorun_ratio;
        goto unlock;
    }
    */

    if (!rqtp) {
        wakeup = p->pcsws.pcsws_job.deadline;
        goto unlock;
    }

  //-  printk(KERN_EMERG "[Wait interval] Setting1\n");

    wakeup = timespec_to_ns(rqtp);

    /*
     * If the task wants to wake up _before_ its absolute deadline
     * we must be sure that reusing its (actual) runtime and deadline
     * at that time _would_ overcome its bandwidth limitation, so
     * that we know it will be given new parameters.
     *
     * If this is not true, we postpone the wake-up time up to the right
     * instant. This involves a division (to calculate the reverse of the
     * task's bandwidth), but it is worth to notice that it is quite
     * unlikely that we get into here very often.
     */
    if (dl_time_before(wakeup, pcsws_se->pcsws_job.deadline) &&
        pcsws_check_bandwidth(pcsws_se, wakeup)) {
        u64 ibw = (u64)pcsws_se->pcsws_job.charge * pcsws_se->pcsws_ded.deadline;

        ibw = div_u64(ibw, pcsws_se->pcsws_ded.budget);

        wakeup = pcsws_se->pcsws_job.deadline - ibw;
      //-  printk(KERN_EMERG "[Wait interval] Setting2\n");
    }

unlock:
    task_rq_unlock(rq, p, &flags);
    lrqtp = ns_to_timespec(wakeup);

  //-  printk(KERN_EMERG "[Wait interval]     [Now] %llu     [Wakeup] %llu\n", rq->clock, wakeup);

    /*
     * If the wakeup instant is in the past we skip hrtimer_nanosleep and
     * keep the task runnable by enqueueing it back, considering that the
     * flag pcsws_new has been set to indicate a new instance.
     */
    if (unlikely(dl_time_before(wakeup, rq->clock_task))) {
        enqueue_task_pcsws(rq, p, 0);
        return 0;
    }

    return hrtimer_nanosleep(&lrqtp, rmtp, HRTIMER_MODE_ABS,
                 CLOCK_MONOTONIC);
}

const struct sched_class pcsws_sched_class = {
    .next               = &rt_sched_class,
    .enqueue_task		= enqueue_task_pcsws,
    .dequeue_task		= dequeue_task_pcsws,
    .yield_task         = yield_task_pcsws,

    .check_preempt_curr	= check_preempt_curr_pcsws,

    .pick_next_task		= pick_next_task_pcsws,
    .put_prev_task		= put_prev_task_pcsws,

    #ifdef CONFIG_SMP
    .select_task_rq		= select_task_rq_pcsws,

    /*
    .set_cpus_allowed       = set_cpus_allowed_pcsws,
    .rq_online              = rq_online_pcsws,
    .rq_offline             = rq_offline_pcsws,
    .pre_schedule		= pre_schedule_pcsws,
    .post_schedule		= post_schedule_pcsws,*/
    .task_woken		= task_woken_pcsws,
    #endif

    .set_curr_task		= set_curr_task_pcsws,
    .task_tick		= task_tick_pcsws,
    .task_fork              = task_fork_pcsws,

    .prio_changed       = prio_changed_pcsws,
    .switched_from		= switched_from_pcsws,
    .switched_to		= switched_to_pcsws,
};

/***************************** INIT RQS  ********************************/
/* Initialize the global runqueue @global_rq */
void init_global_rq(struct pcsws_ready *global_rq)
{
    raw_spin_lock_init(&global_rq->lock);
    global_rq->tasks = RB_ROOT;
    global_rq->nr_running = 0;
    /*global_rq->random = 0;*/
}

/* Initialize the local runqueue @rq */
void init_pcsws_rq(struct pcsws_rq *pcsws_rq, struct rq *rq, struct pcsws_ready *global_rq)
{
    pcsws_rq->pcsws_ready = global_rq;
    pcsws_rq->pjobs = RB_ROOT;

    pcsws_rq->nr_running = 0;
    /*
    pcsws_rq->tot_steals = 0;
    pcsws_rq->tot_enqueues = 0;
    pcsws_rq->tot_picks = 0;
    pcsws_rq->tot_pulls = 0;
    pcsws_rq->tot_cs = 0;
    */
    pcsws_rq->earliest_dl = 0;
}
