#include "Data.h"

#include "TH2F.h"
#include "TH2I.h"
#include "TH1I.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TColor.h"
#include "TBox.h"
#include "TPad.h"
#include "TPaletteAxis.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TGraph2D.h"
#include "TCanvas.h"

#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

Data::Data()
{}

Data::Data(const char* filename)
{
    rootFile = 0;
    T_true = 0;
    T_rec = 0;
    T_proj = 0;
    T_proj_data = 0;

    nChannel_u = 2400;
    nChannel_v = 2400;
    nChannel_w = 3456;
    nTime = 2400;

    rec_cluster_id = new vector<int>;
    rec_x  = new vector<vector<double> >;
    rec_y  = new vector<vector<double> >;
    rec_z  = new vector<vector<double> >;
    rec_dQ = new vector<vector<double> >;
    rec_dx = new vector<vector<double> >;
    rec_L  = new vector<vector<double> >;
    rec_u  = new vector<vector<double> >;
    rec_v  = new vector<vector<double> >;
    rec_w  = new vector<vector<double> >;
    rec_t  = new vector<vector<double> >;

    data_cluster_id = new vector<int>;
    data_channel = new vector<vector<int> >;
    data_time_slice = new vector<vector<int> >;
    data_charge = new vector<vector<int> >;

    LoadData(filename);
    Project();
}



void Data::LoadData(const char* filename)
{
    rootFile = TFile::Open(filename);
    if (!rootFile) {
        string msg = "Unable to open "; msg += filename;
        throw runtime_error(msg.c_str());
    }

    T_true = (TTree*)rootFile->Get("T_true");
    T_rec = (TTree*)rootFile->Get("T_rec");
    T_proj = (TTree*)rootFile->Get("T_proj");
    T_proj_data = (TTree*)rootFile->Get("T_proj_data");

    isData = T_true? false : true;
    cout << "loading data? " << isData << endl;

    LoadRec();
    LoadProj();

}

void Data::LoadRec()
{
    T_rec->SetBranchAddress("rec_x", &rec_x);
    T_rec->SetBranchAddress("rec_y", &rec_y);
    T_rec->SetBranchAddress("rec_z", &rec_z);
    T_rec->SetBranchAddress("rec_dQ", &rec_dQ);
    T_rec->SetBranchAddress("rec_dx", &rec_dx);
    T_rec->SetBranchAddress("rec_L", &rec_L);
    T_rec->SetBranchAddress("rec_cluster_id", &rec_cluster_id);
    T_rec->SetBranchAddress("rec_u", &rec_u);
    T_rec->SetBranchAddress("rec_v", &rec_v);
    T_rec->SetBranchAddress("rec_w", &rec_w);
    T_rec->SetBranchAddress("rec_t", &rec_t);

    T_rec->GetEntry(0);
    nCluster = rec_cluster_id->size();
    cout << " rec cluster id: ";
    for (int i=0; i<nCluster; i++) {
        cout << rec_cluster_id->at(i) << " ";
    }
    cout << endl;

}

void Data::LoadProj()
{
    T_proj_data->SetBranchAddress("cluster_id", &data_cluster_id);
    T_proj_data->SetBranchAddress("channel", &data_channel);
    T_proj_data->SetBranchAddress("time_slice", &data_time_slice);
    T_proj_data->SetBranchAddress("charge", &data_charge);

    T_proj_data->GetEntry(0);
    int size = data_cluster_id->size();
    cout << "data cluster id: ";
    for (int i=0; i<size; i++) {
        int id = data_cluster_id->at(i);
        data_cluster_map[id] = i;
        cout << id << " ";
    }
    cout << endl;
}

void Data::DrawDQDX(int bin, TCanvas* c1, int padNo)
{
    TGraph *g = (TGraph*)gROOT->FindObject("g_dqdx");
    if (g) {
        delete g;
    }
    int size = rec_dQ->at(bin).size();
    g = new TGraph(size);

    for (int i=0; i<size; i++) {
        g->SetPoint(i, rec_L->at(bin).at(i),
            rec_dQ->at(bin).at(i)/1000/rec_dx->at(bin).at(i));
    }
    g->SetName("g_dqdx");
    g->SetTitle(TString::Format("cluster %i", rec_cluster_id->at(bin)));
    g->GetXaxis()->SetTitle("Distance [cm]");
    g->GetYaxis()->SetTitle("dQ/dx [1000 e^{-}/cm]");
    c1->cd(padNo);
    g->Draw("ALP");
    c1->GetPad(padNo)->Modified();
    c1->GetPad(padNo)->Update();
}

void Data::DrawProj(int bin, TCanvas* c1, int padNo)
{
    int cluster_id = rec_cluster_id->at(bin);
    int index = data_cluster_map[cluster_id];

    TH2F *h_proj_u = (TH2F*)gROOT->FindObject("h_proj_u");
    TH2F *h_proj_v = (TH2F*)gROOT->FindObject("h_proj_v");
    TH2F *h_proj_w = (TH2F*)gROOT->FindObject("h_proj_w");

    TGraph *g_rec_u = (TGraph*)gROOT->FindObject("g_rec_u");
    TGraph *g_rec_v = (TGraph*)gROOT->FindObject("g_rec_v");
    TGraph *g_rec_w = (TGraph*)gROOT->FindObject("g_rec_w");
    if (h_proj_u) { delete h_proj_u;}
    if (h_proj_v) { delete h_proj_v;}
    if (h_proj_w) { delete h_proj_w;}
    if (g_rec_u) { delete g_rec_u;}
    if (g_rec_v) { delete g_rec_v;}
    if (g_rec_w) { delete g_rec_w;}

    h_proj_u = new TH2F("h_proj_u", "", nChannel_u, 0-0.5, nChannel_u-0.5, nTime, 1-0.5, nTime+1-0.5);
    h_proj_v = new TH2F("h_proj_v", "", nChannel_v, nChannel_u-0.5, nChannel_u+nChannel_v-0.5, nTime, 1-0.5, nTime+1-0.5);
    h_proj_w = new TH2F("h_proj_w", "", nChannel_w, nChannel_u+nChannel_v-0.5, nChannel_u+nChannel_v+nChannel_w-0.5, nTime, 1-0.5, nTime+1-0.5);

    TH2F *h[3] = {h_proj_u, h_proj_v, h_proj_w};
    int size = data_channel->at(index).size();
    TH2F *hc = 0;
    for (int i=0; i<size; i++) {
        int x = data_channel->at(index).at(i);
        int y = data_time_slice->at(index).at(i);
        int z = data_charge->at(index).at(i);
        if (x<nChannel_u) {
            hc = h_proj_u;
        }
        else if (x<nChannel_u+nChannel_v) {
            hc = h_proj_v;
            x -= nChannel_u;
        }
        else {
            hc = h_proj_w;
            x -= (nChannel_u+nChannel_v);
        }
        hc->SetBinContent(x, y, z);
    }

    int size2 = rec_u->at(bin).size();
    g_rec_u = new TGraph(size2);
    g_rec_v = new TGraph(size2);
    g_rec_w = new TGraph(size2);
    TGraph *g[3] = {g_rec_u, g_rec_v, g_rec_w};

    for (int i=0; i<size2; i++) {
        double u = rec_u->at(bin).at(i);
        double v = rec_v->at(bin).at(i);
        double w = rec_w->at(bin).at(i);
        double t = rec_t->at(bin).at(i);
        // double q = rec_dQ->at(bin).at(i);
        g_rec_u->SetPoint(i, u, t);
        g_rec_v->SetPoint(i, v, t);
        g_rec_w->SetPoint(i, w, t);
    }
    for (int i=0; i<3; i++) {
        h[i]->SetTitle(TString::Format("cluster %i", cluster_id));
        h[i]->GetXaxis()->SetTitle("Channel");
        h[i]->GetYaxis()->SetTitle("Time Slice");
        h[i]->GetZaxis()->SetRangeUser(500, 20000);
        int pad = padNo+i;
        c1->cd(pad);
        h[i]->Draw("colz");
        g[i]->SetLineWidth(2);
        g[i]->SetLineColor(6);
        g[i]->Draw("LPsame");
        c1->GetPad(pad)->Modified();
        c1->GetPad(pad)->Update();
    }


}


void Data::Project()
{

    h_proj[0] = new TH2F("h_proj_u", "", nChannel_u, 1-0.5, nChannel_u-0.5, nTime, 1-0.5, nTime-0.5);
    h_proj[1] = new TH2F("h_proj_v", "", nChannel_v, nChannel_u-0.5, nChannel_u+nChannel_v-0.5, nTime, 1-0.5, nTime-0.5);
    h_proj[2] = new TH2F("h_proj_w", "", nChannel_w, nChannel_u+nChannel_v-0.5, nChannel_u+nChannel_v+nChannel_w-0.5, nTime, 1-0.5, nTime-0.5);

    T_proj_data->Project("h_proj_u", "time_slice:channel", "charge/500", "channel<2400");
    T_proj_data->Project("h_proj_v", "time_slice:channel", "charge/500", "channel<2400 && channel<4800");
    T_proj_data->Project("h_proj_w", "time_slice:channel", "charge/500", "channel>4800");

    //creating boxes
    int threshold = 0;
    TBox *box = 0;
    for (int plane=0; plane<3; plane++) {
        for (int i=1; i<=nChannel_u; i++) {
            for (int j=1; j<=nTime; j++) {
                double content = h_proj[plane]->GetBinContent(i, j);
                if (TMath::Abs(content)>threshold) {
                    box = new TBox(
                        h_proj[plane]->GetXaxis()->GetBinLowEdge(i),
                        h_proj[plane]->GetYaxis()->GetBinLowEdge(j),
                        h_proj[plane]->GetXaxis()->GetBinUpEdge(i),
                        h_proj[plane]->GetYaxis()->GetBinUpEdge(j)
                    );
                    box->SetFillColor(kRed);
                    boxes[plane].push_back(box);
                    box_values[plane].push_back(TMath::Abs(content));
                }
            }
        }
        cout << "plane " << plane << ": " << boxes[plane].size() <<  " boxes created. " << endl;
    }

}

void Data::DrawProjection(int plane)
{

    h_proj[plane]->Clear();
    h_proj[plane]->Draw("colz");
    h_proj[plane]->GetZaxis()->SetRangeUser(10, 20);
    gPad->Update();

    TPaletteAxis *palette = (TPaletteAxis*)h_proj[plane]->GetListOfFunctions()->FindObject("palette");
    int size = boxes[plane].size();
    for (int i=0; i<size; i++) {
        boxes[plane][i]->SetFillColor(palette->GetValueColor(box_values[plane][i]));
        boxes[plane][i]->Draw();
    }
}


Data::~Data()
{
    delete rootFile;
}