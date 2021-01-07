// lowercase interface

extern void spikeinit_(int *spm, int *n, int *klu);
extern void spikeinit_nthread_(int *spm,int *n,int *klu,int *nthread);
extern void spikeinit_default_(int *spm);

extern void dspike_tune_(int *spm);
extern void dspike_gbsv_(int *spm,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *B,int *LDB,int *info);
extern void dspike_gbtrf_(int *spm,int *N,int *KL,int *KU,double *A,int *LDA,double *work,int *info);
extern void dspike_gbtrs_(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *work,double *B,int *LDB);
extern void dspike_gbtrsi_(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,double *C,int *LDC,double *A,int *LDA,double *work,double *B,int *LDB);
extern void dspike_gbtrfp_(int *spm,int *N,int *KL,int *KU,double *A,int *LDA,double *work,int *ipiv, int *info);
extern void dspike_gbtrsp_(int *spm,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *work,int *ipiv,double *B,int *LDB);

extern void zspike_tune_(int *spm);
extern void zspike_gbsv_(int *spm,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *B,int *LDB,int *info);
extern void zspike_gbtrf_(int *spm,int *N,int *KL,int *KU,double *A,int *LDA,double *work,int *info);
extern void zspike_gbtrs_(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *work,double *B,int *LDB);
extern void zspike_gbtrsi_(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,double *C,int *LDC,double *A,int *LDA,double *work,double *B,int *LDB);
extern void zspike_gbtrfp_(int *spm,int *N,int *KL,int *KU,double *A,int *LDA,double *work,int *ipiv, int *info);
extern void zspike_gbtrsp_(int *spm,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *work,int *ipiv,double *B,int *LDB);

extern void cspike_tune_(int *spm);
extern void cspike_gbsv_(int *spm,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *B,int *LDB,int *info);
extern void cspike_gbtrf_(int *spm,int *N,int *KL,int *KU,float *A,int *LDA,float *work,int *info);
extern void cspike_gbtrs_(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *work,float *B,int *LDB);
extern void cspike_gbtrsi_(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,float *C,int *LDC,float *A,int *LDA,float *work,float *B,int *LDB);
extern void cspike_gbtrfp_(int *spm,int *N,int *KL,int *KU,float *A,int *LDA,float *work,int *ipiv,int *info);
extern void cspike_gbtrsp_(int *spm,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *work,int *ipiv,float *B,int *LDB);

extern void sspike_tune_(int *spm);
extern void sspike_gbsv_(int *spm,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *B,int *LDB,int *info);
extern void sspike_gbtrf_(int *spm,int *N,int *KL,int *KU,float *A,int *LDA,float *work,int *info);
extern void sspike_gbtrs_(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *work,float *B,int *LDB);
extern void sspike_gbtrsi_(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,float *C,int *LDC,float *A,int *LDA,float *work,float *B,int *LDB);
extern void sspike_gbtrfp_(int *spm,int *N,int *KL,int *KU,float *A,int *LDA,float *work,int *ipiv,int *info);
extern void sspike_gbtrsp_(int *spm,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *work,int *ipiv,float *B,int *LDB);


void spikeinit(int *spm, int *n, int *klu){
     spikeinit_(spm,n,klu);
}
void spikeinit_default(int *spm){
     spikeinit_default_(spm);
}
void spikeinit_nthread(int *spm, int *n, int *klu, int *nthread){
     spikeinit_nthread_(spm, n, klu, nthread);
}
void dspike_tune(int *spm) {
     dspike_tune_(spm);
}
void dspike_gbsv(int *spm,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *B,int *LDB,int *info){
     dspike_gbsv_(spm,N,KL,KU,NRHS,A,LDA,B,LDB,info);
}
void dspike_gbtrf(int *spm,int *N,int *KL,int *KU,double *A,int *LDA,double *work,int *info){
     dspike_gbtrf_(spm,N,KL,KU,A,LDA,work,info);
}
void dspike_gbtrs(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *work,double *B,int *LDB){
     dspike_gbtrs_(spm,TRANS,N,KL,KU,NRHS,A,LDA,work,B,LDB);
}
void dspike_gbtrsi(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,double *C,int *LDC,double *A,int *LDA,double *work,double *B,int *LDB){
     dspike_gbtrsi_(spm,TRANS,N,KL,KU,NRHS,C,LDC,A,LDA,work,B,LDB);
}
void dspike_gbtrfp(int *spm,int *N,int *KL,int *KU,double *A,int *LDA,double *work,int *ipiv,int *info){
     dspike_gbtrfp_(spm,N,KL,KU,A,LDA,work,ipiv,info);
}
void dspike_gbtrsp(int *spm,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *work,int *ipiv,double *B,int *LDB){
     dspike_gbtrsp_(spm,N,KL,KU,NRHS,A,LDA,work,ipiv,B,LDB);
}


void zspike_tune(int *spm) {
     zspike_tune_(spm);
}
void zspike_gbsv(int *spm,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *B,int *LDB,int *info){
     zspike_gbsv_(spm,N,KL,KU,NRHS,A,LDA,B,LDB,info);
}
void zspike_gbtrf(int *spm,int *N,int *KL,int *KU,double *A,int *LDA,double *work,int *info){
     zspike_gbtrf_(spm,N,KL,KU,A,LDA,work,info);
}
void zspike_gbtrs(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *work,double *B,int *LDB){
     zspike_gbtrs_(spm,TRANS,N,KL,KU,NRHS,A,LDA,work,B,LDB);
}
void zspike_gbtrsi(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,double *C,int *LDC,double *A,int *LDA,double *work,double *B,int *LDB){
     zspike_gbtrsi_(spm,TRANS,N,KL,KU,NRHS,C,LDC,A,LDA,work,B,LDB);
}
void zspike_gbtrfp(int *spm,int *N,int *KL,int *KU,double *A,int *LDA,double *work,int *ipiv,int *info){
     zspike_gbtrfp_(spm,N,KL,KU,A,LDA,work,ipiv,info);
}
void zspike_gbtrsp(int *spm,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *work,int *ipiv,double *B,int *LDB){
     zspike_gbtrsp_(spm,N,KL,KU,NRHS,A,LDA,work,ipiv,B,LDB);
}


void cspike_tune(int *spm) {
     cspike_tune_(spm);
}
void cspike_gbsv(int *spm,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *B,int *LDB,int *info){
     cspike_gbsv_(spm,N,KL,KU,NRHS,A,LDA,B,LDB,info);
}
void cspike_gbtrf(int *spm,int *N,int *KL,int *KU,float *A,int *LDA,float *work,int *info){
     cspike_gbtrf_(spm,N,KL,KU,A,LDA,work,info);
}
void cspike_gbtrs(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *work,float *B,int *LDB){
     cspike_gbtrs_(spm,TRANS,N,KL,KU,NRHS,A,LDA,work,B,LDB);
}
void cspike_gbtrsi(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,float *C,int *LDC,float *A,int *LDA,float *work,float *B,int *LDB){
     cspike_gbtrsi_(spm,TRANS,N,KL,KU,NRHS,C,LDC,A,LDA,work,B,LDB);
}
void cspike_gbtrfp(int *spm,int *N,int *KL,int *KU,float *A,int *LDA,float *work,int *ipiv,int *info){
     cspike_gbtrfp_(spm,N,KL,KU,A,LDA,work,ipiv,info);
}
void cspike_gbtrsp(int *spm,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *work,int *ipiv,float *B,int *LDB){
     cspike_gbtrsp_(spm,N,KL,KU,NRHS,A,LDA,work,ipiv,B,LDB);
}


void sspike_tune(int *spm) {
     sspike_tune_(spm);
}
void sspike_gbsv(int *spm,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *B,int *LDB,int *info){
     sspike_gbsv_(spm,N,KL,KU,NRHS,A,LDA,B,LDB,info);
}
void sspike_gbtrf(int *spm,int *N,int *KL,int *KU,float *A,int *LDA,float *work,int *info){
     sspike_gbtrf_(spm,N,KL,KU,A,LDA,work,info);
}
void sspike_gbtrs(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *work,float *B,int *LDB){
     sspike_gbtrs_(spm,TRANS,N,KL,KU,NRHS,A,LDA,work,B,LDB);
}
void sspike_gbtrsi(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,float *C,int *LDC,float *A,int *LDA,float *work,float *B,int *LDB){
     sspike_gbtrsi_(spm,TRANS,N,KL,KU,NRHS,C,LDC,A,LDA,work,B,LDB);
}
void sspike_gbtrfp(int *spm,int *N,int *KL,int *KU,float *A,int *LDA,float *work,int *ipiv,int *info){
     sspike_gbtrfp_(spm,N,KL,KU,A,LDA,work,ipiv,info);
}
void sspike_gbtrsp(int *spm,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *work,int *ipiv,float *B,int *LDB){
     sspike_gbtrsp_(spm,N,KL,KU,NRHS,A,LDA,work,ipiv,B,LDB);
}

// UPPER CASE INTEFACE

extern void DSPIKE_TUNE_(int *spm);
extern void DSPIKE_GBSV_(int *spm,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *B,int *LDB,int *info);
extern void DSPIKE_GBTRF_(int *spm,int *N,int *KL,int *KU,double *A,int *LDA,double *work,int *info);
extern void DSPIKE_GBTRS_(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *work,double *B,int *LDB);
extern void DSPIKE_GBTRSI_(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,double *C,int *LDC,double *A,int *LDA,double *work,double *B,int *LDB);
extern void DSPIKE_GBTRFP_(int *spm,int *N,int *KL,int *KU,double *A,int *LDA,double *work,int *ipiv, int *info);
extern void DSPIKE_GBTRSP_(int *spm,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *work,int *ipiv,double *B,int *LDB);

extern void ZSPIKE_TUNE_(int *spm);
extern void ZSPIKE_GBSV_(int *spm,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *B,int *LDB,int *info);
extern void ZSPIKE_GBTRF_(int *spm,int *N,int *KL,int *KU,double *A,int *LDA,double *work,int *info);
extern void ZSPIKE_GBTRS_(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *work,double *B,int *LDB);
extern void ZSPIKE_GBTRSI_(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,double *C,int *LDC,double *A,int *LDA,double *work,double *B,int *LDB);
extern void ZSPIKE_GBTRFP_(int *spm,int *N,int *KL,int *KU,double *A,int *LDA,double *work,int *ipiv, int *info);
extern void ZSPIKE_GBTRSP_(int *spm,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *work,int *ipiv,double *B,int *LDB);

extern void CSPIKE_TUNE_(int *spm);
extern void CSPIKE_GBSV_(int *spm,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *B,int *LDB,int *info);
extern void CSPIKE_GBTRF_(int *spm,int *N,int *KL,int *KU,float *A,int *LDA,float *work,int *info);
extern void CSPIKE_GBTRS_(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *work,float *B,int *LDB);
extern void CSPIKE_GBTRSI_(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,float *C,int *LDC,float *A,int *LDA,float *work,float *B,int *LDB);
extern void CSPIKE_GBTRFP_(int *spm,int *N,int *KL,int *KU,float *A,int *LDA,float *work,int *ipiv,int *info);
extern void CSPIKE_GBTRSP_(int *spm,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *work,int *ipiv,float *B,int *LDB);

extern void SSPIKE_TUNE_(int *spm);
extern void SSPIKE_GBSV_(int *spm,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *B,int *LDB,int *info);
extern void SSPIKE_GBTRF_(int *spm,int *N,int *KL,int *KU,float *A,int *LDA,float *work,int *info);
extern void SSPIKE_GBTRS_(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *work,float *B,int *LDB);
extern void SSPIKE_GBTRSI_(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,float *C,int *LDC,float *A,int *LDA,float *work,float *B,int *LDB);
extern void SSPIKE_GBTRFP_(int *spm,int *N,int *KL,int *KU,float *A,int *LDA,float *work,int *ipiv,int *info);
extern void SSPIKE_GBTRSP_(int *spm,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *work,int *ipiv,float *B,int *LDB);


void SPIKEINIT(int *spm, int *n, int *klu){
     spikeinit_(spm,n,klu);
}
void SPIKEINIT_DEFAULT(int *spm){
     spikeinit_default_(spm);
}
void SPIKEINIT_NTHREAD(int *spm, int *n, int *klu, int *nthread){
     spikeinit_nthread_(spm, n, klu, nthread);
}
void DSPIKE_TUNE(int *spm) {
     dspike_tune_(spm);
}
void DSPIKE_GBSV(int *spm,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *B,int *LDB,int *info){
     dspike_gbsv_(spm,N,KL,KU,NRHS,A,LDA,B,LDB,info);
}
void DSPIKE_GBTRF(int *spm,int *N,int *KL,int *KU,double *A,int *LDA,double *work,int *info){
     dspike_gbtrf_(spm,N,KL,KU,A,LDA,work,info);
}
void DSPIKE_GBTRS(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *work,double *B,int *LDB){
     dspike_gbtrs_(spm,TRANS,N,KL,KU,NRHS,A,LDA,work,B,LDB);
}
void DSPIKE_GBTRSI(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,double *C,int *LDC,double *A,int *LDA,double *work,double *B,int *LDB){
     dspike_gbtrsi_(spm,TRANS,N,KL,KU,NRHS,C,LDC,A,LDA,work,B,LDB);
}
void DSPIKE_GBTRFP(int *spm,int *N,int *KL,int *KU,double *A,int *LDA,double *work,int *ipiv,int *info){
     dspike_gbtrfp_(spm,N,KL,KU,A,LDA,work,ipiv,info);
}
void DSPIKE_GBTRSP(int *spm,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *work,int *ipiv,double *B,int *LDB){
     dspike_gbtrsp_(spm,N,KL,KU,NRHS,A,LDA,work,ipiv,B,LDB);
}


void ZSPIKE_TUNE(int *spm) {
     zspike_tune_(spm);
}
void ZSPIKE_GBSV(int *spm,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *B,int *LDB,int *info){
     zspike_gbsv_(spm,N,KL,KU,NRHS,A,LDA,B,LDB,info);
}
void ZSPIKE_GBTRF(int *spm,int *N,int *KL,int *KU,double *A,int *LDA,double *work,int *info){
     zspike_gbtrf_(spm,N,KL,KU,A,LDA,work,info);
}
void ZSPIKE_GBTRS(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *work,double *B,int *LDB){
     zspike_gbtrs_(spm,TRANS,N,KL,KU,NRHS,A,LDA,work,B,LDB);
}
void ZSPIKE_GBTRSI(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,double *C,int *LDC,double *A,int *LDA,double *work,double *B,int *LDB){
     zspike_gbtrsi_(spm,TRANS,N,KL,KU,NRHS,C,LDC,A,LDA,work,B,LDB);
}
void ZSPIKE_GBTRFP(int *spm,int *N,int *KL,int *KU,double *A,int *LDA,double *work,int *ipiv,int *info){
     zspike_gbtrfp_(spm,N,KL,KU,A,LDA,work,ipiv,info);
}
void ZSPIKE_GBTRSP(int *spm,int *N,int *KL,int *KU,int *NRHS,double *A,int *LDA,double *work,int *ipiv,double *B,int *LDB){
     zspike_gbtrsp_(spm,N,KL,KU,NRHS,A,LDA,work,ipiv,B,LDB);
}


void CSPIKE_TUNE(int *spm) {
     cspike_tune_(spm);
}
void CSPIKE_GBSV(int *spm,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *B,int *LDB,int *info){
     cspike_gbsv_(spm,N,KL,KU,NRHS,A,LDA,B,LDB,info);
}
void CSPIKE_GBTRF(int *spm,int *N,int *KL,int *KU,float *A,int *LDA,float *work,int *info){
     cspike_gbtrf_(spm,N,KL,KU,A,LDA,work,info);
}
void CSPIKE_GBTRS(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *work,float *B,int *LDB){
     cspike_gbtrs_(spm,TRANS,N,KL,KU,NRHS,A,LDA,work,B,LDB);
}
void CSPIKE_GBTRSI(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,float *C,int *LDC,float *A,int *LDA,float *work,float *B,int *LDB){
     cspike_gbtrsi_(spm,TRANS,N,KL,KU,NRHS,C,LDC,A,LDA,work,B,LDB);
}
void CSPIKE_GBTRFP(int *spm,int *N,int *KL,int *KU,float *A,int *LDA,float *work,int *ipiv,int *info){
     cspike_gbtrfp_(spm,N,KL,KU,A,LDA,work,ipiv,info);
}
void CSPIKE_GBTRSP(int *spm,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *work,int *ipiv,float *B,int *LDB){
     cspike_gbtrsp_(spm,N,KL,KU,NRHS,A,LDA,work,ipiv,B,LDB);
}


void SSPIKE_TUNE(int *spm) {
     sspike_tune_(spm);
}
void SSPIKE_GBSV(int *spm,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *B,int *LDB,int *info){
     sspike_gbsv_(spm,N,KL,KU,NRHS,A,LDA,B,LDB,info);
}
void SSPIKE_GBTRF(int *spm,int *N,int *KL,int *KU,float *A,int *LDA,float *work,int *info){
     sspike_gbtrf_(spm,N,KL,KU,A,LDA,work,info);
}
void SSPIKE_GBTRS(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *work,float *B,int *LDB){
     sspike_gbtrs_(spm,TRANS,N,KL,KU,NRHS,A,LDA,work,B,LDB);
}
void SSPIKE_GBTRSI(int *spm,char *TRANS,int *N,int *KL,int *KU,int *NRHS,float *C,int *LDC,float *A,int *LDA,float *work,float *B,int *LDB){
     sspike_gbtrsi_(spm,TRANS,N,KL,KU,NRHS,C,LDC,A,LDA,work,B,LDB);
}
void SSPIKE_GBTRFP(int *spm,int *N,int *KL,int *KU,float *A,int *LDA,float *work,int *ipiv,int *info){
     sspike_gbtrfp_(spm,N,KL,KU,A,LDA,work,ipiv,info);
}
void SSPIKE_GBTRSP(int *spm,int *N,int *KL,int *KU,int *NRHS,float *A,int *LDA,float *work,int *ipiv,float *B,int *LDB){
     sspike_gbtrsp_(spm,N,KL,KU,NRHS,A,LDA,work,ipiv,B,LDB);
}

