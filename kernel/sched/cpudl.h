#ifndef _LINUX_CPUDL_H
#define _LINUX_CPUDL_H

#include <linux/sched.h>

#define IDX_INVALID     -1

struct array_item {
    u64 dl;
    int cpu;
};

struct cpudl {
    raw_spinlock_t lock;
    int size;
    int cpu_to_idx[NR_CPUS];
    struct array_item elements[NR_CPUS];
    cpumask_var_t free_cpus;
};


#ifdef CONFIG_SMP
int find_idle_cpu_pcsws(struct cpudl *cp);
int find_latest_cpu_pcsws(struct cpudl *cp);
int find_random_stealable_cpu_pcsws(struct cpudl *cp, int cpu);
int find_earliest_stealable_cpu_pcsws(struct cpudl *cp);
void cpudl_set(struct cpudl *cp, int cpu, u64 dl, int is_valid);
int cpudl_init(struct cpudl *cp);
void cpudl_cleanup(struct cpudl *cp);
inline int cpudl_maximum(struct cpudl *cp);

#else
#define cpudl_set(cp, cpu, dl) do { } while (0)
#define cpudl_init() do { } while (0)
#endif /* CONFIG_SMP */

#endif /* _LINUX_CPUDL_H */
