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
#include "TMarker.h"

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
    c1 = 0;
    pad_dqdx = 1;
    pad_proj = 4;
    currentCluster = 0;
    doDrawBadCh = false;

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

    for (int i=0; i<3; i++) {
        currentPoint[i] = new TMarker(0, 0, 24);
        currentPoint[i]->SetMarkerColor(6);
        currentPoint[i]->SetMarkerSize(4);
    }

    LoadData(filename);
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
    T_bad_ch = (TTree*)rootFile->Get("T_bad_ch");

    isData = T_true? false : true;
    cout << "loading data? " << isData << endl;

    LoadRec();
    LoadProj();
    LoadBadCh();

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

void Data::LoadBadCh()
{
    if (!T_bad_ch) {
        cout << "no bad channles in the file ..." << endl;
        return;
    }
    int chid, start, end;
    TLine *line = 0;
    T_bad_ch->SetBranchAddress("chid", &chid);
    T_bad_ch->SetBranchAddress("start_time", &start);
    T_bad_ch->SetBranchAddress("end_time", &end);
    int nEntries = T_bad_ch->GetEntries();
    double rebin = 4;
    for (int i=0; i<nEntries; i++) {
        T_bad_ch->GetEntry(i);
        bad_id.push_back(chid);
        bad_start.push_back(start);
        bad_end.push_back(end);
        line = new TLine(chid, start/rebin, chid, end/rebin);
        line->SetLineColorAlpha(kGray, 0.5);
        bad_lines.push_back(line);
    }
}

void Data::DrawDQDX()
{
    TGraph *g = (TGraph*)gROOT->FindObject("g_dqdx");
    if (g) {
        delete g;
    }
    int size = rec_dQ->at(currentCluster).size();
    g = new TGraph(size);

    for (int i=0; i<size; i++) {
        g->SetPoint(i, rec_L->at(currentCluster).at(i),
            rec_dQ->at(currentCluster).at(i)/1000/rec_dx->at(currentCluster).at(i));
    }
    g->SetName("g_dqdx");
    g->SetTitle(TString::Format("cluster %i", rec_cluster_id->at(currentCluster)));
    g->GetXaxis()->SetTitle("Distance [cm]");
    g->GetYaxis()->SetTitle("dQ/dx [1000 e^{-}/cm]");
    c1->cd(pad_dqdx);
    g->Draw("ALP");
    c1->GetPad(pad_dqdx)->Modified();
    c1->GetPad(pad_dqdx)->Update();
}

void Data::DrawProj()
{
    int cluster_id = rec_cluster_id->at(currentCluster);
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

    h_proj_u = new TH2F("h_proj_u", "", nChannel_u, 0, nChannel_u, nTime, 0, nTime);
    h_proj_v = new TH2F("h_proj_v", "", nChannel_v, nChannel_u, nChannel_u+nChannel_v, nTime, 0, nTime);
    h_proj_w = new TH2F("h_proj_w", "", nChannel_w, nChannel_u+nChannel_v, nChannel_u+nChannel_v+nChannel_w, nTime, 0, nTime);

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
        hc->SetBinContent(x+1, y+1, z);
    }

    int size2 = rec_u->at(currentCluster).size();
    g_rec_u = new TGraph(size2);
    g_rec_v = new TGraph(size2);
    g_rec_w = new TGraph(size2);
    TGraph *g[3] = {g_rec_u, g_rec_v, g_rec_w};

    for (int i=0; i<size2; i++) {
        double u = rec_u->at(currentCluster).at(i);
        double v = rec_v->at(currentCluster).at(i);
        double w = rec_w->at(currentCluster).at(i);
        double t = rec_t->at(currentCluster).at(i);
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
        int pad = pad_proj+i;
        c1->cd(pad);
        h[i]->Draw("colz");
        g[i]->SetLineWidth(2);
        g[i]->SetLineColor(6);
        g[i]->Draw("LPsame");
        c1->GetPad(pad)->Modified();
        c1->GetPad(pad)->Update();
    }

    DrawBadCh();
}

void Data::ZoomProj(int pointIndex, int zoomBin)
{
    int u =  rec_u->at(currentCluster).at(pointIndex);
    int v =  rec_v->at(currentCluster).at(pointIndex);
    int w =  rec_w->at(currentCluster).at(pointIndex);
    int t =  rec_t->at(currentCluster).at(pointIndex);

    cout << " index: " << pointIndex;
    cout << " u: " << u;
    cout << " v: " << v;
    cout << " w: " << w;
    cout << " t: " << t;
    cout << endl;

    TH2F *h = (TH2F*)gROOT->FindObject("h_proj_u");
    h->GetXaxis()->SetRangeUser(u-zoomBin, u+zoomBin);
    h->GetYaxis()->SetRangeUser(t-zoomBin, t+zoomBin);
    c1->GetPad(pad_proj)->Modified();
    c1->GetPad(pad_proj)->Update();

    h = (TH2F*)gROOT->FindObject("h_proj_v");
    h->GetXaxis()->SetRangeUser(v-zoomBin, v+zoomBin);
    h->GetYaxis()->SetRangeUser(t-zoomBin, t+zoomBin);
    c1->GetPad(pad_proj+1)->Modified();
    c1->GetPad(pad_proj+1)->Update();

    h = (TH2F*)gROOT->FindObject("h_proj_w");
    h->GetXaxis()->SetRangeUser(w-zoomBin, w+zoomBin);
    h->GetYaxis()->SetRangeUser(t-zoomBin, t+zoomBin);
    c1->GetPad(pad_proj+2)->Modified();
    c1->GetPad(pad_proj+2)->Update();
}

void Data::DrawPoint(int pointIndex)
{
    double x[3] = {0};
    x[0] =  rec_u->at(currentCluster).at(pointIndex);
    x[1] =  rec_v->at(currentCluster).at(pointIndex);
    x[2] =  rec_w->at(currentCluster).at(pointIndex);
    double t =  rec_t->at(currentCluster).at(pointIndex);

    for (int i=0; i<3; i++) {
        c1->cd(pad_proj+i);
        currentPoint[i]->SetX(x[i]);
        currentPoint[i]->SetY(t);
        currentPoint[i]->Draw();
        c1->GetPad(pad_proj+i)->Modified();
        c1->GetPad(pad_proj+i)->Update();
    }

}

void Data::DrawBadCh()
{
    int size = bad_lines.size();
    for (int pad=pad_proj; pad<pad_proj+3; pad++) {
        c1->cd(pad);
        TVirtualPad *p = c1->GetPad(pad);
        for (int i=0; i<size; i++) {
            if ( (pad==pad_proj && bad_id[i]<2400)
                || (pad==pad_proj+1 && bad_id[i]>=2400 && bad_id[i]<4800)
                || (pad==pad_proj+2 && bad_id[i]>=4800)
            ) {
                if (doDrawBadCh) {
                    bad_lines[i]->Draw();
                }
                else {
                   p->GetListOfPrimitives()->Remove(bad_lines[i]);
                }

            }
        }
        p->Modified();
        p->Update();
    }

}




// void Data::Project()
// {

//     h_proj[0] = new TH2F("h_proj_u", "", nChannel_u, 1-0.5, nChannel_u-0.5, nTime, 1-0.5, nTime-0.5);
//     h_proj[1] = new TH2F("h_proj_v", "", nChannel_v, nChannel_u-0.5, nChannel_u+nChannel_v-0.5, nTime, 1-0.5, nTime-0.5);
//     h_proj[2] = new TH2F("h_proj_w", "", nChannel_w, nChannel_u+nChannel_v-0.5, nChannel_u+nChannel_v+nChannel_w-0.5, nTime, 1-0.5, nTime-0.5);

//     T_proj_data->Project("h_proj_u", "time_slice:channel", "charge/500", "channel<2400");
//     T_proj_data->Project("h_proj_v", "time_slice:channel", "charge/500", "channel<2400 && channel<4800");
//     T_proj_data->Project("h_proj_w", "time_slice:channel", "charge/500", "channel>4800");

//     //creating boxes
//     int threshold = 0;
//     TBox *box = 0;
//     for (int plane=0; plane<3; plane++) {
//         for (int i=1; i<=nChannel_u; i++) {
//             for (int j=1; j<=nTime; j++) {
//                 double content = h_proj[plane]->GetBinContent(i, j);
//                 if (TMath::Abs(content)>threshold) {
//                     box = new TBox(
//                         h_proj[plane]->GetXaxis()->GetBinLowEdge(i),
//                         h_proj[plane]->GetYaxis()->GetBinLowEdge(j),
//                         h_proj[plane]->GetXaxis()->GetBinUpEdge(i),
//                         h_proj[plane]->GetYaxis()->GetBinUpEdge(j)
//                     );
//                     box->SetFillColor(kRed);
//                     boxes[plane].push_back(box);
//                     box_values[plane].push_back(TMath::Abs(content));
//                 }
//             }
//         }
//         cout << "plane " << plane << ": " << boxes[plane].size() <<  " boxes created. " << endl;
//     }

// }

// void Data::DrawProjection(int plane)
// {

//     h_proj[plane]->Clear();
//     h_proj[plane]->Draw("colz");
//     h_proj[plane]->GetZaxis()->SetRangeUser(10, 20);
//     gPad->Update();

//     TPaletteAxis *palette = (TPaletteAxis*)h_proj[plane]->GetListOfFunctions()->FindObject("palette");
//     int size = boxes[plane].size();
//     for (int i=0; i<size; i++) {
//         boxes[plane][i]->SetFillColor(palette->GetValueColor(box_values[plane][i]));
//         boxes[plane][i]->Draw();
//     }
// }


Data::~Data()
{
    delete rootFile;
}