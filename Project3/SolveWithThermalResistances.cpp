vector<double> calculateMultiLayerTempDistribution(double missionDuration, double dt,
    const vector<int>& nodesPerLayer,
    const vector<double>& layerDx_m,
    const vector<MaterialProperties>& layerMaterials) {

    // Step 1: Build the mesh
    vector<double> x_positions;
    vector<double> k;      // thermal conductivity at each node
    vector<double> rho_cp; // rho * cp at each node

    for (size_t i = 0; i < nodesPerLayer.size(); ++i) {
        for (int j = 0; j < nodesPerLayer[i]; ++j) {
            double xpos = (j * layerDx_m[i]);
            if (!x_positions.empty()) xpos += x_positions.back();
            x_positions.push_back(xpos);

            k.push_back(layerMaterials[i].thermalConductivity);
            rho_cp.push_back(layerMaterials[i].density * layerMaterials[i].specificHeatCapacity);
        }
        if (!x_positions.empty()) x_positions.pop_back(); // remove duplicate interface node
        if (!k.empty()) k.pop_back();
        if (!rho_cp.empty()) rho_cp.pop_back();
    }

    int N = x_positions.size();
    vector<double> T(N, 300.0); // Initial temperature
    T[0] = 900.0;               // Boundary condition: fixed 900K at outer surface

    // Step 2: Setup tridiagonal matrix using THERMAL RESISTANCES
    vector<double> a(N-1, 0.0); // sub-diagonal
    vector<double> b(N, 0.0);   // main diagonal
    vector<double> c(N-1, 0.0); // super-diagonal

    double dx = (x_positions[1] - x_positions[0]); // assume locally constant dx

    // Compute thermal resistances at faces
    vector<double> R_face(N-1, 0.0);
    for (int i = 0; i < N-1; ++i) {
        R_face[i] = dx / (2.0 * k[i]) + dx / (2.0 * k[i+1]);
    }

    // Assemble matrix using R_face
    for (int i = 1; i < N-1; ++i) {
        double coeff_left = dt / (rho_cp[i] * R_face[i-1] * dx);
        double coeff_right = dt / (rho_cp[i] * R_face[i] * dx);

        a[i-1] = -coeff_left;
        b[i] = 1.0 + coeff_left + coeff_right;
        c[i] = -coeff_right;
    }

    // Boundary conditions
    b[0] = 1.0;     
    c[0] = 0.0;

    double last_coeff = dt / (rho_cp[N-1] * R_face[N-2] * dx);
    a[N-2] = -2.0 * last_coeff;
    b[N-1] = 1.0 + 2.0 * last_coeff; // insulated boundary

    // Step 3: Time stepping loop
    for (double time = 0; time < missionDuration; time += dt) {
        vector<double> d = T;
        d[0] = 900.0; 
        T = thomas_algorithm(a, b, c, d);
    }

    return T;
}
