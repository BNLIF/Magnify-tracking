#ifndef DATA_H
#define DATA_H

#include <vector>
#include <map>

class TFile;
class TTree;
class TH2F;
class TBox;
class TCanvas;

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


    TH2F *h_proj[3];
    std::vector<TBox*> boxes[3];
    vector<float> box_values[3];


    bool isData;
    int nChannel_u;
    int nChannel_v;
    int nChannel_w;
    int nTime;

    void Project();
    void LoadRec();
    void LoadProj();
    void DrawProjection(int plane);
    void DrawDQDX(int bin, TCanvas* c1, int padNo=1);
    void DrawProj(int bin, TCanvas* c1, int padNo=4);

private:
    void LoadData(const char* filename);


};

#endif
