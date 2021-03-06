#include "GuiController.h"
#include "MainWindow.h"
#include "ViewWindow.h"
#include "ControlWindow.h"
#include "Data.h"

#include "TApplication.h"
#include "TSystem.h"
#include "TExec.h"
#include "TROOT.h"
#include "TMath.h"
#include "TGFileDialog.h"

#include "TGMenu.h"
#include "TGNumberEntry.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TH2I.h"
#include "TH1I.h"
#include "TBox.h"
#include "TLine.h"
#include "TColor.h"
#include "TLatex.h"

#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


GuiController::GuiController(const TGWindow *p, int w, int h, const char* fn, int sign)
{
    mw = new MainWindow(p, w, h);
    vw = mw->fViewWindow;
    cw = mw->fControlWindow;
    currentCluster = 0;
    currentPointIndex = 0;

    TString filename;
    if (!fn) {
        filename = OpenDialog();
    }
    else {
        filename = fn;
    }
    data = new Data(filename.Data(), sign);
    data->c1 = vw->can;
    SetCurrentCluster(0);
    data->DrawNewCluster();


    // mw->SetWindowName(TString::Format("Magnify: run %i, sub-run %i, event %i",
    //     data->runNo, data->subRunNo, data->eventNo));

    // for (int i=0; i<6; i++) {
    //     vw->can->cd(i+1);
    //     data->wfs.at(i)->Draw2D();
    // }
    // for (int i=0; i<3; i++) {
    //     vw->can->cd(i+7);
    //     int chanNo = data->wfs.at(i)->firstChannel;
    //     data->wfs.at(i)->Draw1D(chanNo);
    //     TH1F *h = data->wfs.at(i+3)->Draw1D(chanNo, "same"); // draw calib
    //     h->SetLineColor(kRed);
    //     hCurrent[i] = h;
    // }



    InitConnections();

}

GuiController::~GuiController()
{
    // gApplication->Terminate(0);
}

void GuiController::InitConnections()
{
    mw->fMenuFile->Connect("Activated(int)", "GuiController", this, "HandleMenu(int)");

    cw->clusterEntry->SetLimitValues(0, data->nCluster-1);
    cw->clusterEntry->Connect("ValueSet(Long_t)", "GuiController", this, "ClusterChanged(int)");

    cw->clusterIdEntry->SetNumber(data->rec_cluster_id->at(currentCluster));
    cw->clusterIdEntry->Connect("ValueSet(Long_t)", "GuiController", this, "ClusterIdChanged(int)");

    cw->pointIndexEntry->Connect("ValueSet(Long_t)", "GuiController", this, "ZoomChanged()");

    cw->zoomEntry->Connect("ValueSet(Long_t)", "GuiController", this, "ZoomChanged()");
    cw->badChanelButton->Connect("Clicked()", "GuiController", this, "ToggleBadChannel()");
    cw->drawTrackButton->Connect("Clicked()", "GuiController", this, "ToggleDrawTrack()");
    cw->allClusterButton->Connect("Clicked()", "GuiController", this, "ToggleAllCluster()");

    cw->unZoomButton->Connect("Clicked()", "GuiController", this, "UnZoom()");
    cw->rangeZoomButton->Connect("Clicked()", "GuiController", this, "RangeZoom()");

    vw->can->Connect(
        "ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",
        "GuiController",
        this,
        "ProcessCanvasEvent(Int_t,Int_t,Int_t,TObject*)"
    );
}

void GuiController::ToggleBadChannel()
{
    data->doDrawBadCh = cw->badChanelButton->IsDown();
    data->DrawBadCh();
}

void GuiController::ToggleDrawTrack()
{
    data->doDrawTrack = cw->drawTrackButton->IsDown();
    data->DrawProj();
    ZoomChanged();
}

void GuiController::ToggleAllCluster()
{
    if (cw->allClusterButton->IsDown()) {
        data->DrawProjAll(
            cw->minEntry[0]->GetNumber(), cw->maxEntry[0]->GetNumber(),
            cw->minEntry[1]->GetNumber(), cw->maxEntry[1]->GetNumber(),
            cw->minEntry[2]->GetNumber(), cw->maxEntry[2]->GetNumber(),
            cw->minEntry[3]->GetNumber(), cw->maxEntry[3]->GetNumber()
        );
        // ZoomChanged();
    }
    else {
        data->DrawProj();
        ZoomChanged();

    }
;
}


void GuiController::SetCurrentCluster(int newCluster)
{
    cout << "new cluster: " << newCluster << endl;
    currentCluster = newCluster;
    data->currentCluster = currentCluster;

    int size = data->rec_dQ->at(currentCluster).size();
    cw->pointIndexEntry->SetLimitValues(0, size-1);
    cw->pointIndexEntry->SetNumber(0);
}

void GuiController::ClusterChanged(int i)
{
    int newCluster = cw->clusterEntry->GetNumber();
    if (newCluster == currentCluster) return;
    SetCurrentCluster(newCluster);

    cw->clusterIdEntry->SetNumber(data->rec_cluster_id->at(newCluster));

    data->DrawNewCluster();
    if (cw->keepRangeButton->IsDown()) {
        RangeZoom();
    }

}

void GuiController::ClusterIdChanged(int i)
{
    int newId = cw->clusterIdEntry->GetNumber();
    int newCluster = data->rec_cluster_map[newId];
    // cout << newCluster << endl;

    cw->clusterEntry->SetNumber(newCluster);

    if (newCluster == currentCluster) return;
    SetCurrentCluster(newCluster);

    data->DrawNewCluster();
    if (cw->keepRangeButton->IsDown()) {
        RangeZoom();
    }
}


void GuiController::UnZoom()
{
    data->ZoomProj(0, -1);
    SetRangeEntries();
}

void GuiController::RangeZoom()
{
    data->ZoomProj(0, -1,
        cw->minEntry[0]->GetNumber(), cw->maxEntry[0]->GetNumber(),
        cw->minEntry[1]->GetNumber(), cw->maxEntry[1]->GetNumber(),
        cw->minEntry[2]->GetNumber(), cw->maxEntry[2]->GetNumber(),
        cw->minEntry[3]->GetNumber(), cw->maxEntry[3]->GetNumber()
    );
}

void GuiController::ZoomChanged()
{
    currentPointIndex = cw->pointIndexEntry->GetNumber();
    int zoomBin = 10*cw->zoomEntry->GetNumber();
    data->ZoomProj(currentPointIndex, zoomBin);
    data->DrawPoint(currentPointIndex);
    SetRangeEntries();
}

void GuiController::ProcessCanvasEvent(Int_t ev, Int_t x, Int_t y, TObject *selected)
{
    if (ev == 11) { // clicked
        if (!(selected->IsA() == TGraph::Class()
            || selected->IsA() == TH2F::Class()
            // || selected->IsA() == TLine::Class()
        )) return;
        TVirtualPad* pad = vw->can->GetClickSelectedPad();
        int padNo = pad->GetNumber();
        double xx = pad->AbsPixeltoX(x);
        double yy = pad->AbsPixeltoY(y);
        int ci = data->FindClusterIndex(xx, yy);
        int pi = data->FindPointIndex(xx, yy);
        cout << "pad " << padNo << ": (" << xx << ", " << yy << ")" 
            << "; cluster index: " << ci
            << "; point index: " << pi
            << endl;
        TString text = TString::Format("(%.0f, %.0f): point index %d", xx, yy, pi);
        data->infoText->SetText(data->infoText->GetX(), data->infoText->GetY(), text.Data());
        TVirtualPad *infoPad = vw->can->GetPad(data->pad_dqdx+1);
        infoPad->cd();
        infoPad->Modified();
        infoPad->Update();

        if (selected->IsA() == TGraph::Class() && padNo <=3) { // first row tgraph clicked
            TGraph *g = (TGraph*)gROOT->FindObject("g_dqdx");
            currentPointIndex = TMath::BinarySearch(g->GetN(), g->GetX(), xx);
    	    if (currentPointIndex!=g->GetN()-1){
    	      if (fabs(xx-*(g->GetX()+currentPointIndex)) >
    		  fabs(xx-*(g->GetX()+currentPointIndex+1)))
    		  currentPointIndex += 1;
    	    }
            int zoomBin = 10*cw->zoomEntry->GetNumber();

	    // std::cout << zoomBin << " " << cw->zoomEntry->GetNumber() << std::endl;
            cw->pointIndexEntry->SetNumber(currentPointIndex);

            data->ZoomProj(currentPointIndex, zoomBin);
            data->DrawPoint(currentPointIndex);
            SetRangeEntries();
        }
        else if (selected->IsA() == TH2F::Class()) { // 2nd row th2f clicked

        }
    }
}

void GuiController::SetRangeEntries()
{    
    cw->minEntry[0]->SetNumber(vw->can->GetPad(data->pad_proj)->GetUymin());
    cw->maxEntry[0]->SetNumber(vw->can->GetPad(data->pad_proj)->GetUymax());
    cw->minEntry[1]->SetNumber(vw->can->GetPad(data->pad_proj)->GetUxmin());
    cw->maxEntry[1]->SetNumber(vw->can->GetPad(data->pad_proj)->GetUxmax());
    cw->minEntry[2]->SetNumber(vw->can->GetPad(data->pad_proj+1)->GetUxmin());
    cw->maxEntry[2]->SetNumber(vw->can->GetPad(data->pad_proj+1)->GetUxmax());
    cw->minEntry[3]->SetNumber(vw->can->GetPad(data->pad_proj+2)->GetUxmin());
    cw->maxEntry[3]->SetNumber(vw->can->GetPad(data->pad_proj+2)->GetUxmax());
}

void GuiController::HandleMenu(int id)
{
    // const char *filetypes[] = {"ROOT files", "*.root", 0, 0};
    switch (id) {
        case M_FILE_EXIT:
            gApplication->Terminate(0);
            break;
    }
}

TString GuiController::OpenDialog()
{
    const char *filetypes[] = {"ROOT files", "*.root", 0, 0};
    TString currentDir(gSystem->WorkingDirectory());
    static TString dir(".");
    TGFileInfo fi;
    fi.fFileTypes = filetypes;
    fi.fIniDir    = StrDup(dir);
    new TGFileDialog(gClient->GetRoot(), mw, kFDOpen, &fi);
    dir = fi.fIniDir;
    gSystem->cd(currentDir.Data());

    if (fi.fFilename) {
        // UnZoom();
        cout << "open file: " << fi.fFilename << endl;
        return fi.fFilename;
    }
    else {
        gApplication->Terminate(0);
    }
    return "";

}
