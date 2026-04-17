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
class TH3F;
class TPolyLine3D;
class TPolyMarker3D;
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

    std::vector<int> *rec_cluster_id;
    std::vector<std::vector<double> >* rec_x;
    std::vector<std::vector<double> >* rec_y;
    std::vector<std::vector<double> >* rec_z;
    std::vector<std::vector<double> >* rec_dQ;
    std::vector<std::vector<double> >* rec_dx;
    std::vector<std::vector<double> >* rec_L;
    std::vector<std::vector<double> >* rec_u;
    std::vector<std::vector<double> >* rec_v;
    std::vector<std::vector<double> >* rec_w;
    std::vector<std::vector<double> >* rec_t;
    std::vector<std::vector<double> >* reduced_chi2;
    std::vector<std::vector<double> >* com_dis;
    std::vector<std::vector<double> >* com_dtheta;
    std::vector<std::vector<double> >* true_dQ;
    std::vector<std::vector<int> >* flag_vertex;
    std::vector<std::vector<int> >* sub_cluster_id;
    int nCluster;
    std::map<int, int> rec_cluster_map; // cluster id => index
    std::vector<int> sub_id; // list of sub cluster IDs.
    std::vector<int> sub_start_index;
    std::vector<int> sub_end_index;

    std::vector<int> *data_cluster_id;
    std::vector<std::vector<int> >* data_channel;
    std::vector<std::vector<int> >* data_time_slice;
    std::vector<std::vector<int> >* data_charge;
    std::vector<std::vector<int> >* data_charge_pred;
    std::vector<std::vector<int> >* data_charge_err;
    std::map<int, int> data_cluster_map; // cluster id => index

    std::vector<int> bad_id;
    std::vector<int> bad_start;
    std::vector<int> bad_end;

    std::vector<std::vector<double> >* true_x;
    std::vector<std::vector<double> >* true_y;
    std::vector<std::vector<double> >* true_z;

    std::vector<double>* stat_beg_dis;
    std::vector<double>* stat_end_dis;



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
    std::vector<TLine*> bad_lines;
    std::vector<TGraph*> g_subclusters_u;
    std::vector<TGraph*> g_subclusters_v;
    std::vector<TGraph*> g_subclusters_w;
    std::vector<TPolyLine3D*>   g_subclusters_3d_lines;
    std::vector<TPolyMarker3D*> g_subclusters_3d_marks;
    TH3F          *h_3d_frame;
    TPolyMarker3D *g_3d_points;
    TPolyLine3D   *gt_3d_line;
    TPolyMarker3D *g_3d_start;
    TPolyMarker3D *currentPoint3d;
    TLatex *infoText;

    // Persistent projection / prediction histograms. Allocated once via
    // EnsureProjHistos() and Reset() between cluster changes.
    TH2F *h_proj_u;
    TH2F *h_proj_v;
    TH2F *h_proj_w;
    TH2F *h_pred_u;
    TH2F *h_pred_v;
    TH2F *h_pred_w;
    TH2F *h_proj_u_all;
    TH2F *h_proj_v_all;
    TH2F *h_proj_w_all;
    TH2F *h_pred_u_all;
    TH2F *h_pred_v_all;
    TH2F *h_pred_w_all;

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
    void EnsureProjHistos();
    void EnsureProjAllHistos();
    static double applySign(double residual, int sign);


};

#endif
