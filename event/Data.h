#ifndef DATA_H
#define DATA_H

#include <vector>
#include <map>

class TFile;
class TTree;
class TH2F;
class TBox;
class TCanvas;
class TMarker;
class TLine;
class TGraph;
class TGraph2D;
class TLatex;

class Data {
public:
    Data();
    Data(const char* filename, int sign=0);

    virtual ~Data();

    TFile *rootFile;
    int sign;

    TTree *T_true;
    TTree *T_rec;
    TTree *T_proj;
    TTree *T_proj_data;
    TTree *T_bad_ch;

    vector<int> *rec_cluster_id;
    vector<vector<double> >* rec_x;
    vector<vector<double> >* rec_y;
    vector<vector<double> >* rec_z;
    vector<vector<double> >* rec_dQ;
    vector<vector<double> >* rec_dx;
    vector<vector<double> >* rec_L;
    vector<vector<double> >* rec_u;
    vector<vector<double> >* rec_v;
    vector<vector<double> >* rec_w;
    vector<vector<double> >* rec_t;
    vector<vector<double> >* reduced_chi2;
    vector<vector<double> >* com_dis;
    vector<vector<double> >* com_dtheta;
    vector<vector<double> >* true_dQ;
    vector<vector<int> >* flag_vertex;
    vector<vector<int> >* sub_cluster_id;
    int nCluster;
    std::map<int, int> rec_cluster_map; // cluster id => index
    vector<int> sub_id; // list of sub cluster IDs.
    vector<int> sub_start_index;
    vector<int> sub_end_index; 

    vector<int> *data_cluster_id;
    vector<vector<int> >* data_channel;
    vector<vector<int> >* data_time_slice;
    vector<vector<int> >* data_charge;
    vector<vector<int> >* data_charge_pred;
    vector<vector<int> >* data_charge_err;
    std::map<int, int> data_cluster_map; // cluster id => index

    vector<int> bad_id;
    vector<int> bad_start;
    vector<int> bad_end;

    vector<vector<double> >* true_x;
    vector<vector<double> >* true_y;
    vector<vector<double> >* true_z;

    vector<double>* stat_beg_dis;
    vector<double>* stat_end_dis;



    // TH2F *h_proj[3];
    // std::vector<TBox*> boxes[3];
    // vector<float> box_values[3];


    bool isData;
    int nChannel_u;
    int nChannel_v;
    int nChannel_w;
    int nTime;

    TCanvas *c1;
    int pad_proj;
    int pad_dqdx;
    int pad_mc_compare;
    int pad_3d;
    int pad_pred;
    int currentCluster;
    bool doDrawBadCh;
    bool doDrawTrack;
    TMarker *currentPoint[3];
    TMarker *dqdxPoint;
    vector<TLine*> bad_lines;
    vector<TGraph*> g_subclusters_u;
    vector<TGraph*> g_subclusters_v;
    vector<TGraph*> g_subclusters_w;
    vector<TGraph2D*> g_subclusters_3d;
    TLatex *infoText;

    void LoadRec();
    void LoadProj();
    void LoadBadCh();
    void LoadTruth();

    void DrawDQDX();
    void DrawProj();    
    void DrawProjAll(int t0=-2, int t1=-1, int u0=-1, int u1=-1, int v0=-1, int v1=-1, int w0=-1, int w1=-1);
    void ZoomProj(int pointIndex, int zoomBin, 
        int t0=-2, int t1=-1, int u0=-1, int u1=-1, int v0=-1, int v1=-1, int w0=-1, int w1=-1);
    void DrawSubclusters();
    void DrawPoint(int pointIndex);
    void DrawBadCh();
    void Draw3D();
    void DrawMCCompare();

    void DrawNewCluster();
    void Clear();

    int FindClusterIndex(double x, double y);
    int FindPointIndex(double x, double y);


private:
    void LoadData(const char* filename);


};

#endif
