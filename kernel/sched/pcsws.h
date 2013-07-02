#ifndef _LINUX_PCSWS_H
#define _LINUX_PCSWS_H

enum pcsws_event{
    ENQUEUE=0,
    DEQUEUE,
    SWITCH_TO,
    SWITCH_AWAY,
    SERVER_EXPIRED,
    SERVER_EXHAUSTED,
};

void init_timer_pcsws(struct sched_pcsws_entity *pcsws_se);
inline u64 get_next_activation_period(struct sched_pcsws_entity *pcsws_se);
inline int deadline_missed_pcsws(struct rq *rq, struct sched_pcsws_entity *pcsws_se);
inline void update_task_pcsws(struct rq *rq, struct sched_pcsws_entity *pcsws_se);
void start_timer_pcsws(struct rq *rq, struct sched_pcsws_entity *pcsws_se, u64 diff);
long wait_interval_pcsws(struct task_struct *p, struct timespec *rqtp,
              struct timespec *rmtp);

void init_global_rq(struct pcsws_ready *global_rq);
void init_pcsws_rq(struct pcsws_rq *pcsws_rq, struct rq *rq, struct pcsws_ready *global_rq);
void init_pcsws_job(struct rq *rq, struct task_struct *p);

#endif
