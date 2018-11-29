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

class Data {
public:
    Data();
    Data(const char* filename);

    virtual ~Data();

    TFile *rootFile;

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
    int nCluster;

    vector<int> *data_cluster_id;
    vector<vector<int> >* data_channel;
    vector<vector<int> >* data_time_slice;
    vector<vector<int> >* data_charge;
    std::map<int, int> data_cluster_map; // cluster id => index

    vector<int> bad_id;
    vector<int> bad_start;
    vector<int> bad_end;



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
    int pad_3d;
    int currentCluster;
    bool doDrawBadCh;
    TMarker *currentPoint[3];
    vector<TLine*> bad_lines;

    void LoadRec();
    void LoadProj();
    void LoadBadCh();

    void DrawDQDX();
    void DrawProj();
    void ZoomProj(int pointIndex, int zoomBin);
    void DrawPoint(int pointIndex);
    void DrawBadCh();
    void Draw3D();

    void DrawNewCluster();


private:
    void LoadData(const char* filename);


};

#endif
