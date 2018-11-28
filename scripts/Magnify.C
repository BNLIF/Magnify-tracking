void Magnify(const char* filename="../data/track_com_data.root");
void Magnify(const char* filename=0)
{
    // Data *data = new Data(filename);

    GuiController *gc = new GuiController(
        gClient->GetRoot(),
        1600,
        900,
        filename
    );


}