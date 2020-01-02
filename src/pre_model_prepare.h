
/* flag of media_parameterization */
#define LOC 0
#define GRI 1
#define VOL 2
#define ORT 3
#define TTI 4

#define PREPARE_ISO 0
#define PREPARE_ORT 1
#define PREPARE_TTI 2

int model_prepare(int nbmodels, char *materialfile, char *interfacesfile,
        int effective_para_method, int save_model, int *prepared_media,
        float *xvec1, float *zvec1, float *xvec2, float *zvec2, int nx, int nz, float dx, float dz,
        float *c11_1, float *c13_1, float *c15_1, float *c33_1, float *c35_1, float *c55_1,
        float *c11_2, float *c13_2, float *c15_2, float *c33_2, float *c35_2, float *c55_2,
        float *B01, float *B10,
        float *lam2mu00, float *lam2mu11, float *lam00, float *lam11, float *mu00, float *mu11,
        float *rho01, float *rho10, float *rho_x, float *rho_z);

int write_model_iso(int save_model, float *lam2mu, float *lam, float *muxz,
        float *rho_x, float *rho_z, int nx_all, int nz_all);

int write_model_TTI(int save_model,
        float *c11_1, float *c13_1, float *c33_1,
        float *c15_1, float *c35_1, float *c55_1,
        float *c11_2, float *c13_2, float *c33_2,
        float *c15_2, float *c35_2, float *c55_2,
        float *B01, float *B10, int nx_all, int nz_all);

//int boundary_expansion(float *u, int nx, int nz,
//        int nghost_x1, int nghost_x2, int nghost_z1, int nghost_z2);
