#define R0 ((i_loc-1)   * dr)
#define R1 ((i_loc+1-1) * dr)
#define R2 ((i_loc+2-1) * dr)
#define R3 ((i_loc+3-1) * dr)
#define DR3 (dr*dr*dr)

#define ALPHA_0(r)   ((R3 - r)*(R2 - r)*(R1 - r)/(6*DR3))
#define ALPHA_1(r)   ((R3 - r)*(R2 - r)*(r - R0)/(2*DR3))
#define ALPHA_2(r)   ((R3 - r)*(r - R1)*(r - R0)/(2*DR3))
#define ALPHA_3(r)   ((r - R2)*(r - R1)*(r - R0)/(6*DR3))

#define D_ALPHA_0(r) (( -(R2 - (r))*(R1 - (r)) - (R3 - (r))*(R1 - (r)) - (R3 - (r))*(R2 - (r)) ) \
		      /(6*DR3))
#define D_ALPHA_1(r) (( -(R2 - (r))*((r) - R0) - (R3 - (r))*((r) - R0) + (R3 - (r))*(R2 - (r)) ) \
		      /(2*DR3))
#define D_ALPHA_2(r) (( -((r) - R1)*((r) - R0) + (R3 - (r))*((r) - R0) + (R3 - (r))*((r) - R1) ) \
		      /(2*DR3))
#define D_ALPHA_3(r) (( +((r) - R1)*((r) - R0) + ((r) - R2)*((r) - R0) + ((r) - R2)*((r) - R1) ) \
		      /(6*DR3))

#define D2_ALPHA_0(r) ((-6.*(r) + 2.*R1 + 2.*R2 + 2.*R3 )/(6*DR3))
#define D2_ALPHA_1(r) (( 6.*(r) - 2.*R0 - 2.*R2 - 2.*R3 )/(2*DR3))
#define D2_ALPHA_2(r) ((-6.*(r) + 2.*R0 + 2.*R1 + 2.*R3 )/(2*DR3))
#define D2_ALPHA_3(r) (( 6.*(r) - 2.*R1 - 2.*R0 - 2.*R2 )/(6*DR3))

#define S0 ((j_loc-1)  * ds)
#define S1 ((j_loc+1-1)* ds)
#define S2 ((j_loc+2-1) * ds)
#define S3 ((j_loc+3-1) * ds)
#define DS3 (ds*ds*ds)

#define BETA_0(s)   ((S3 - s)*(S2 - s)*(S1 - s)/(6*DS3))
#define BETA_1(s)   ((S3 - s)*(S2 - s)*(s - S0)/(2*DS3))
#define BETA_2(s)   ((S3 - s)*(s - S1)*(s - S0)/(2*DS3))
#define BETA_3(s)   ((s - S2)*(s - S1)*(s - S0)/(6*DS3))

#define D_BETA_0(s) (( -(S2 - (s))*(S1 - (s)) - (S3 - (s))*(S1 - (s)) - (S3 - (s))*(S2 - (s)) ) \
		      /(6*DS3))
#define D_BETA_1(s) (( -(S2 - (s))*((s) - S0) - (S3 - (s))*((s) - S0) + (S3 - (s))*(S2 - (s)) ) \
		      /(2*DS3))
#define D_BETA_2(s) (( -((s) - S1)*((s) - S0) + (S3 - (s))*((s) - S0) + (S3 - (s))*((s) - S1) ) \
		      /(2*DS3))
#define D_BETA_3(s) (( +((s) - S1)*((s) - S0) + ((s) - S2)*((s) - S0) + ((s) - S2)*((s) - S1) ) \
		      /(6*DS3))

#define D2_BETA_0(s) ((-6.*(s) + 2.*S1 + 2.*S2 + 2.*S3 )/(6*DS3))
#define D2_BETA_1(s) (( 6.*(s) - 2.*S0 - 2.*S2 - 2.*S3 )/(2*DS3))
#define D2_BETA_2(s) ((-6.*(s) + 2.*S0 + 2.*S1 + 2.*S3 )/(2*DS3))
#define D2_BETA_3(s) (( 6.*(s) - 2.*S1 - 2.*S0 - 2.*S2 )/(6*DS3))

#define CUB(r,s,x) (ALPHA_0(r) * BETA_0(s) * x(i_m    ,j_m    ) +  \
		    ALPHA_1(r) * BETA_0(s) * x(i_loc+1,j_m    ) +  \
		    ALPHA_2(r) * BETA_0(s) * x(i_loc+2,j_m    ) +  \
		    ALPHA_3(r) * BETA_0(s) * x(i_p    ,j_m    ) +  \
		    ALPHA_0(r) * BETA_1(s) * x(i_m    ,j_loc+1) +  \
		    ALPHA_1(r) * BETA_1(s) * x(i_loc+1,j_loc+1) +  \
		    ALPHA_2(r) * BETA_1(s) * x(i_loc+2,j_loc+1) +  \
		    ALPHA_3(r) * BETA_1(s) * x(i_p    ,j_loc+1) +  \
		    ALPHA_0(r) * BETA_2(s) * x(i_m    ,j_loc+2) +  \
		    ALPHA_1(r) * BETA_2(s) * x(i_loc+1,j_loc+2) +  \
		    ALPHA_2(r) * BETA_2(s) * x(i_loc+2,j_loc+2) +  \
		    ALPHA_3(r) * BETA_2(s) * x(i_p    ,j_loc+2) +  \
		    ALPHA_0(r) * BETA_3(s) * x(i_m    ,j_p    ) +  \
		    ALPHA_1(r) * BETA_3(s) * x(i_loc+1,j_p    ) +  \
		    ALPHA_2(r) * BETA_3(s) * x(i_loc+2,j_p    ) +  \
		    ALPHA_3(r) * BETA_3(s) * x(i_p    ,j_p    ) )

#define CUB_R(r,s,x) (D_ALPHA_0(r) * BETA_0(s) * x(i_m    ,j_m    ) +  \
		      D_ALPHA_1(r) * BETA_0(s) * x(i_loc+1,j_m    ) +  \
		      D_ALPHA_2(r) * BETA_0(s) * x(i_loc+2,j_m    ) +  \
		      D_ALPHA_3(r) * BETA_0(s) * x(i_p    ,j_m    ) +  \
		      D_ALPHA_0(r) * BETA_1(s) * x(i_m    ,j_loc+1) +  \
		      D_ALPHA_1(r) * BETA_1(s) * x(i_loc+1,j_loc+1) +  \
		      D_ALPHA_2(r) * BETA_1(s) * x(i_loc+2,j_loc+1) +  \
		      D_ALPHA_3(r) * BETA_1(s) * x(i_p    ,j_loc+1) +  \
		      D_ALPHA_0(r) * BETA_2(s) * x(i_m    ,j_loc+2) +  \
		      D_ALPHA_1(r) * BETA_2(s) * x(i_loc+1,j_loc+2) +  \
		      D_ALPHA_2(r) * BETA_2(s) * x(i_loc+2,j_loc+2) +  \
		      D_ALPHA_3(r) * BETA_2(s) * x(i_p    ,j_loc+2) +  \
		      D_ALPHA_0(r) * BETA_3(s) * x(i_m    ,j_p    ) +  \
		      D_ALPHA_1(r) * BETA_3(s) * x(i_loc+1,j_p    ) +  \
		      D_ALPHA_2(r) * BETA_3(s) * x(i_loc+2,j_p    ) +  \
		      D_ALPHA_3(r) * BETA_3(s) * x(i_p    ,j_p    ) )

#define CUB_S(r,s,x) (ALPHA_0(r) * D_BETA_0(s) * x(i_m    ,j_m    ) +  \
		      ALPHA_1(r) * D_BETA_0(s) * x(i_loc+1,j_m    ) +  \
		      ALPHA_2(r) * D_BETA_0(s) * x(i_loc+2,j_m    ) +  \
		      ALPHA_3(r) * D_BETA_0(s) * x(i_p    ,j_m    ) +  \
		      ALPHA_0(r) * D_BETA_1(s) * x(i_m    ,j_loc+1) +  \
		      ALPHA_1(r) * D_BETA_1(s) * x(i_loc+1,j_loc+1) +  \
		      ALPHA_2(r) * D_BETA_1(s) * x(i_loc+2,j_loc+1) +  \
		      ALPHA_3(r) * D_BETA_1(s) * x(i_p    ,j_loc+1) +  \
		      ALPHA_0(r) * D_BETA_2(s) * x(i_m    ,j_loc+2) +  \
		      ALPHA_1(r) * D_BETA_2(s) * x(i_loc+1,j_loc+2) +  \
		      ALPHA_2(r) * D_BETA_2(s) * x(i_loc+2,j_loc+2) +  \
		      ALPHA_3(r) * D_BETA_2(s) * x(i_p    ,j_loc+2) +  \
		      ALPHA_0(r) * D_BETA_3(s) * x(i_m    ,j_p    ) +  \
		      ALPHA_1(r) * D_BETA_3(s) * x(i_loc+1,j_p    ) +  \
		      ALPHA_2(r) * D_BETA_3(s) * x(i_loc+2,j_p    ) +  \
		      ALPHA_3(r) * D_BETA_3(s) * x(i_p    ,j_p    ) )

#define CUB_RR(r,s,x) (D2_ALPHA_0(r) * BETA_0(s) * x(i_m    ,j_m    ) +  \
		       D2_ALPHA_1(r) * BETA_0(s) * x(i_loc+1,j_m    ) +  \
		       D2_ALPHA_2(r) * BETA_0(s) * x(i_loc+2,j_m    ) +  \
		       D2_ALPHA_3(r) * BETA_0(s) * x(i_p    ,j_m    ) +  \
		       D2_ALPHA_0(r) * BETA_1(s) * x(i_m    ,j_loc+1) +  \
		       D2_ALPHA_1(r) * BETA_1(s) * x(i_loc+1,j_loc+1) +  \
		       D2_ALPHA_2(r) * BETA_1(s) * x(i_loc+2,j_loc+1) +  \
		       D2_ALPHA_3(r) * BETA_1(s) * x(i_p    ,j_loc+1) +  \
		       D2_ALPHA_0(r) * BETA_2(s) * x(i_m    ,j_loc+2) +  \
		       D2_ALPHA_1(r) * BETA_2(s) * x(i_loc+1,j_loc+2) +  \
		       D2_ALPHA_2(r) * BETA_2(s) * x(i_loc+2,j_loc+2) +  \
		       D2_ALPHA_3(r) * BETA_2(s) * x(i_p    ,j_loc+2) +  \
		       D2_ALPHA_0(r) * BETA_3(s) * x(i_m    ,j_p    ) +  \
		       D2_ALPHA_1(r) * BETA_3(s) * x(i_loc+1,j_p    ) +  \
		       D2_ALPHA_2(r) * BETA_3(s) * x(i_loc+2,j_p    ) +  \
		       D2_ALPHA_3(r) * BETA_3(s) * x(i_p    ,j_p    ) )

#define CUB_RS(r,s,x) (D_ALPHA_0(r) * D_BETA_0(s) * x(i_m    ,j_m    ) +  \
		       D_ALPHA_1(r) * D_BETA_0(s) * x(i_loc+1,j_m    ) +  \
		       D_ALPHA_2(r) * D_BETA_0(s) * x(i_loc+2,j_m    ) +  \
		       D_ALPHA_3(r) * D_BETA_0(s) * x(i_p    ,j_m    ) +  \
		       D_ALPHA_0(r) * D_BETA_1(s) * x(i_m    ,j_loc+1) +  \
		       D_ALPHA_1(r) * D_BETA_1(s) * x(i_loc+1,j_loc+1) +  \
		       D_ALPHA_2(r) * D_BETA_1(s) * x(i_loc+2,j_loc+1) +  \
		       D_ALPHA_3(r) * D_BETA_1(s) * x(i_p    ,j_loc+1) +  \
		       D_ALPHA_0(r) * D_BETA_2(s) * x(i_m    ,j_loc+2) +  \
		       D_ALPHA_1(r) * D_BETA_2(s) * x(i_loc+1,j_loc+2) +  \
		       D_ALPHA_2(r) * D_BETA_2(s) * x(i_loc+2,j_loc+2) +  \
		       D_ALPHA_3(r) * D_BETA_2(s) * x(i_p    ,j_loc+2) +  \
		       D_ALPHA_0(r) * D_BETA_3(s) * x(i_m    ,j_p    ) +  \
		       D_ALPHA_1(r) * D_BETA_3(s) * x(i_loc+1,j_p    ) +  \
		       D_ALPHA_2(r) * D_BETA_3(s) * x(i_loc+2,j_p    ) +  \
		       D_ALPHA_3(r) * D_BETA_3(s) * x(i_p    ,j_p    ) )

#define CUB_SS(r,s,x) (ALPHA_0(r) * D2_BETA_0(s) * x(i_m    ,j_m    ) +  \
		       ALPHA_1(r) * D2_BETA_0(s) * x(i_loc+1,j_m    ) +  \
		       ALPHA_2(r) * D2_BETA_0(s) * x(i_loc+2,j_m    ) +  \
		       ALPHA_3(r) * D2_BETA_0(s) * x(i_p    ,j_m    ) +  \
		       ALPHA_0(r) * D2_BETA_1(s) * x(i_m    ,j_loc+1) +  \
		       ALPHA_1(r) * D2_BETA_1(s) * x(i_loc+1,j_loc+1) +  \
		       ALPHA_2(r) * D2_BETA_1(s) * x(i_loc+2,j_loc+1) +  \
		       ALPHA_3(r) * D2_BETA_1(s) * x(i_p    ,j_loc+1) +  \
		       ALPHA_0(r) * D2_BETA_2(s) * x(i_m    ,j_loc+2) +  \
		       ALPHA_1(r) * D2_BETA_2(s) * x(i_loc+1,j_loc+2) +  \
		       ALPHA_2(r) * D2_BETA_2(s) * x(i_loc+2,j_loc+2) +  \
		       ALPHA_3(r) * D2_BETA_2(s) * x(i_p    ,j_loc+2) +  \
		       ALPHA_0(r) * D2_BETA_3(s) * x(i_m    ,j_p    ) +  \
		       ALPHA_1(r) * D2_BETA_3(s) * x(i_loc+1,j_p    ) +  \
		       ALPHA_2(r) * D2_BETA_3(s) * x(i_loc+2,j_p    ) +  \
		       ALPHA_3(r) * D2_BETA_3(s) * x(i_p    ,j_p    ) )
