#ifndef CHP_PROC_H
#define CHP_PROC_H

typedef struct chp_proc_ chp_proc;
struct chp_proc_ {
    int rank;
    int group_size;
};

void chp_proc_init(chp_proc *P);
void chp_proc_fini(chp_proc *P);

#endif // CHP_PROC_H
