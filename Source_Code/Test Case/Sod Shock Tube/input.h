#define initial_dx 0.0025
#define Co 0.5
#define max_dt initial_dx*Co
#define ref_iter 10
#define end_time 0.1
#define max_iter end_time/(max_dt*ref_iter)
#define min_div 8
#define max_refinement 4
#define Lx 1
#define Ly 0.06
#define dimX Lx/(initial_dx*min_div)
#define dimY Ly/(initial_dx*min_div)
#define dimZ dimY

#define r_l 1.0
#define P_l 0.0
#define r_r 0.125
#define P_r 0.1

#define gam 1.4
