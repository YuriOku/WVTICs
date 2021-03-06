float betamodel ( const float r );
float GalaxyCluster_Density ( const int ipart , const double bias );

void GalaxyCluster_Velocity ( const int ipart, double out[3] );
float GalaxyCluster_U ( const int ipart );

extern struct HaloProperties {
    float Rho0;
    float Beta;
    float Rcore;
    float R_Sample;
} Halo;
