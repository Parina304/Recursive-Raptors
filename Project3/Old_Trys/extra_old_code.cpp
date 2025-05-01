//Generates dx values per layer based on number of nodes 
vector<double> generateDxList(const vector<double>& layerThicknesses_cm, const vector<int>& nodes_per_layer) {
    vector<double> dx_per_layer_m;

    for (size_t i = 0; i < layerThicknesses_cm.size(); ++i) {
        double dx = layerThicknesses_cm[i] / 100.0 / (nodes_per_layer[i] - 1); // -1 because we include start and end node & div by 100 for cm to m  
        dx_per_layer_m.push_back(dx);
    }

    return dx_per_layer_m;
}

//Resolves dx values per layer such that the nodes are equally spaced 
pair<vector<double>, vector<int>> resolveDxList(const vector<double>& layerThicknesses_cm, const vector<double>& dx_guess_per_layer_cm) {
    
    vector<double> resolved_dx_per_layer_m;
    vector<int> nodes_per_layer;

    for (size_t i = 0; i < layerThicknesses_cm.size(); ++i) {
        double thickness = layerThicknesses_cm[i] / 100.0; // cm → meters
        double dx_guess = dx_guess_per_layer_cm[i] / 100.0; // cm → meters

        // Calculate approximate number of intervals
        int intervals = (int) round(thickness / dx_guess);
        if (intervals < 1) intervals = 1; // Avoid zero division
            
        double dx_exact = thickness / intervals;

        resolved_dx_per_layer_m.push_back(dx_exact);
        nodes_per_layer.push_back(intervals + 1); // Nodes = intervals + 1
    }

    return {resolved_dx_per_layer_m, nodes_per_layer};
}