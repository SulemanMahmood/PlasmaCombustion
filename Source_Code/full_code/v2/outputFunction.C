// this is a member function of Cell Chare

void Cell::writeOutput() {
    char filename[100];
    sprintf(filename, "CellOutput_%d_%d_%d_%d.dat", thisIndex.x, thisIndex.y, thisIndex.z, iter); // iter may need to be change to meet the real name
    FILE* OutFile = fopen(filename, "w");
    fprintf(OutFile, "VARIABLES = \"VX\" \"VY\" \"VZ\" \"Density\" \"Temp\" \"CH4\" \"O2\" \"OH\" \"CH3\" \"H\" \"O\" \"CO2\" \"H2O\"\n");
    fprintf(OutFile, "ZONE I=%d, J=%d, K=%d, DATAPACKING=POINT\n", ndiv, ndiv, ndiv);
    
    for (int k = 0; k < ndiv; k++) {
        for (int j = 0; j < ndiv; j++) {
            for (int i = 0; i < ndiv; i++) {
                fprintf(OutFile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", val_new[i][j][k].u, val_new[i][j][k].v, val_new[i][j][k].w, val_new[i][j][k].r, Tg[i][j][k], val_new[i][j][k].Y[10], val_new[i][j][k].Y[3], val_new[i][j][k].Y[4], val_new[i][j][k].Y[9], val_new[i][j][k].Y[1], val_new[i][j][k].Y[2], val_new[i][j][k].Y[12], val_new[i][j][k].Y[5]);
            }
        }
    } // end triple loop
}