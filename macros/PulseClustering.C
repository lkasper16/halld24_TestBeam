
TString treeName = "gem_trackers_hits";
int minimum_NPulse = 50;

struct Point {
    float x, y, weight;
    bool visited;
    int clusterId;
    Point(float x, float y, float weight) : x(x), y(y), weight(weight), visited(false), clusterId(-1) {}
};

float distance(const Point &a, const Point &b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

void expandCluster(vector<Point> &points_x, Point &point, int clusterId, float eps, int minPts) {
    vector<Point*> seeds;
    for (auto &p : points_x) {
        if (distance(point, p) <= eps) {
            seeds.push_back(&p);
        }
    }

    if (seeds.size() < minPts) {
        point.clusterId = 0; // Mark as noise
        return;
    }

    for (auto &seed : seeds) {
        seed->clusterId = clusterId;
    }

    seeds.erase(remove(seeds.begin(), seeds.end(), &point), seeds.end());

    while (!seeds.empty()) {
        Point *current = seeds.back();
        seeds.pop_back();

        if (!current->visited) {
            current->visited = true;

            vector<Point*> result;
            for (auto &p : points_x) {
                if (distance(*current, p) <= eps) {
                    result.push_back(&p);
                }
            }

            if (result.size() >= minPts) {
                for (auto &res : result) {
                    if (res->clusterId == -1 || res->clusterId == 0) {
                        if (res->clusterId == -1) {
                            seeds.push_back(res);
                        }
                        res->clusterId = clusterId;
                    }
                }
            }
        }
    }
}

void dbscan(vector<Point> &points_x, float eps, int minPts) {
    int clusterId = 1;
    for (auto &point : points_x) {
        if (!point.visited) {
            point.visited = true;
            expandCluster(points_x, point, clusterId, eps, minPts);
            if (point.clusterId == clusterId) {
                clusterId++;
            }
        }
    }
}

void PulseClustering(TString fileName = "/work/halld2/home/nseptian/halld24_TestBeam/run/test_4349_thr150_wRadiators_clustering.root"){
    gStyle->SetOptStat(0);
    TFile *f = TFile::Open(fileName);
    TTree *t = (TTree*)f->Get(treeName);

    cout << "Number of entries: " << t->GetEntries() << endl;

    vector<int> *v_pulsecut_xpos = 0;
    vector<float> *v_pulsecut_xpos_time = 0;
    vector<float> *v_pulsecut_xpos_amp = 0;
    vector<int> *v_pulsecut_ypos = 0;
    vector<float> *v_pulsecut_ypos_time = 0;
    vector<float> *v_pulsecut_ypos_amp = 0;

    t->SetBranchAddress("v_pulsecut_xpos", &v_pulsecut_xpos);
    t->SetBranchAddress("v_pulsecut_xpos_time", &v_pulsecut_xpos_time);
    t->SetBranchAddress("v_pulsecut_xpos_amp", &v_pulsecut_xpos_amp);
    t->SetBranchAddress("v_pulsecut_ypos", &v_pulsecut_ypos);
    t->SetBranchAddress("v_pulsecut_ypos_time", &v_pulsecut_ypos_time);
    t->SetBranchAddress("v_pulsecut_ypos_amp", &v_pulsecut_ypos_amp);

    TH2F *h_time_xpos_pulsecut = new TH2F("h_time_xpos_pulsecut", "", 200, 0, 200, 121, -0.5, 120.5);
    TH2F *h_time_ypos_pulsecut = new TH2F("h_time_ypos_pulsecut", "", 200, 0, 200, 528, -0.5, 527.5);

    TH2F *h_cluster_xpos = new TH2F("h_cluster_xpos", "", 200, 0, 200, 121, -0.5, 120.5);
    TH2F *h_cluster_ypos = new TH2F("h_cluster_ypos", "", 200, 0, 200, 528, -0.5, 527.5);
    
    TCanvas *c = new TCanvas("c", "c", 1600, 800);
    c->Divide(2, 1);

    TString pdfFileName = "clustering_output_10events.pdf";
    c->Print(pdfFileName + "[");

    Int_t filled_entry = 0;
    for (int i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);

        if ((v_pulsecut_xpos->size() < minimum_NPulse) || (v_pulsecut_ypos->size() < minimum_NPulse)) continue;

        h_time_xpos_pulsecut->Reset();
        h_time_ypos_pulsecut->Reset();

        for (size_t j = 0; j < v_pulsecut_xpos->size(); ++j) {
            h_time_xpos_pulsecut->Fill(v_pulsecut_xpos_time->at(j), v_pulsecut_xpos->at(j), v_pulsecut_xpos_amp->at(j));
        }

        for (size_t j = 0; j < v_pulsecut_ypos->size(); ++j) {
            h_time_ypos_pulsecut->Fill(v_pulsecut_ypos_time->at(j), v_pulsecut_ypos->at(j), v_pulsecut_ypos_amp->at(j));
        }

        vector<Point> points_x;
        for (int x = 1; x <= h_time_xpos_pulsecut->GetNbinsX(); ++x) {
            for (int y = 1; y <= h_time_xpos_pulsecut->GetNbinsY(); ++y) {
                float binContent = h_time_xpos_pulsecut->GetBinContent(x, y);
                if (binContent > 0) {
                    points_x.emplace_back(h_time_xpos_pulsecut->GetXaxis()->GetBinCenter(x), h_time_xpos_pulsecut->GetYaxis()->GetBinCenter(y), binContent);
                }
            }
        }

        vector<Point> points_y;
        for (int x = 1; x <= h_time_ypos_pulsecut->GetNbinsX(); ++x) {
            for (int y = 1; y <= h_time_ypos_pulsecut->GetNbinsY(); ++y) {
                float binContent = h_time_ypos_pulsecut->GetBinContent(x, y);
                if (binContent > 0) {
                    points_y.emplace_back(h_time_ypos_pulsecut->GetXaxis()->GetBinCenter(x), h_time_ypos_pulsecut->GetYaxis()->GetBinCenter(y), binContent);
                }
            }
        }

        float eps = 5.0; // Set the epsilon value for DBSCAN (maximum distance between two points to be considered as in the same neighborhood)
        int minPts = 2; // Set the minimum number of points_x for DBSCAN
        dbscan(points_x, eps, minPts);
        dbscan(points_y, eps, minPts);

        h_cluster_xpos->Reset();
        for (const auto &point : points_x) {
            if (point.clusterId > 0) { // Ignore noise points_x
                cout << "(x) Cluster ID: " << point.clusterId << " time: " << point.x << " x: " << point.y << endl;
            }
        }

        map<int, pair<float, float>> clusterSums_x;
        map<int, int> clusterCounts_x;

        for (const auto &point : points_x) {
            if (point.clusterId > 0) { // Ignore noise points_x
                clusterSums_x[point.clusterId].first += point.x*point.weight;
                clusterSums_x[point.clusterId].second += point.y*point.weight;
                clusterCounts_x[point.clusterId]+=point.weight;
            }
        }

        for (const auto &cluster : clusterSums_x) {
            int clusterId = cluster.first;
            float centroidX = cluster.second.first / clusterCounts_x[clusterId];
            float centroidY = cluster.second.second / clusterCounts_x[clusterId];
            h_cluster_xpos->Fill(centroidX, centroidY);
        }

        map<int, pair<float, float>> clusterSums_y;
        map<int, int> clusterCounts_y;

        for (const auto &point : points_y) {
            if (point.clusterId > 0) { // Ignore noise points_x
                clusterSums_y[point.clusterId].first += point.x*point.weight;
                clusterSums_y[point.clusterId].second += point.y*point.weight;
                clusterCounts_y[point.clusterId]+=point.weight;
            }
        }

        h_cluster_ypos->Reset();
        for (const auto &point : points_y) {
            if (point.clusterId > 0) { // Ignore noise points_x
                cout << "(y) Cluster ID: " << point.clusterId << " time: " << point.x << " y: " << point.y << endl;
            }
        }

        for (const auto &cluster : clusterSums_y) {
            int clusterId = cluster.first;
            float centroidX = cluster.second.first / clusterCounts_y[clusterId];
            float centroidY = cluster.second.second / clusterCounts_y[clusterId];
            h_cluster_ypos->Fill(centroidX, centroidY);
        }

        c->cd(1);
        h_time_xpos_pulsecut->Draw("COLZ");
        h_time_xpos_pulsecut->GetXaxis()->SetTitle("Time (sample)");
        h_time_xpos_pulsecut->GetYaxis()->SetTitle("X channel (mm)");
        h_time_xpos_pulsecut->SetMinimum(150);
        
        h_cluster_xpos->Draw("SAME");
        h_cluster_xpos->SetMarkerStyle(kCircle);
        h_cluster_xpos->SetMarkerSize(2.0);

        c->cd(2);
        h_time_ypos_pulsecut->Draw("COLZ");
        h_time_ypos_pulsecut->GetXaxis()->SetTitle("Time (sample)");
        h_time_ypos_pulsecut->GetYaxis()->SetTitle("Y channel (mm)");
        h_time_ypos_pulsecut->SetMinimum(150);

        h_cluster_ypos->Draw("SAME");
        h_cluster_ypos->SetMarkerStyle(kCircle);
        h_cluster_ypos->SetMarkerSize(2.0);

        c->Print(pdfFileName);
        filled_entry++;
        if (filled_entry > 10) break;
    }

    c->Print(pdfFileName + "]");

}
