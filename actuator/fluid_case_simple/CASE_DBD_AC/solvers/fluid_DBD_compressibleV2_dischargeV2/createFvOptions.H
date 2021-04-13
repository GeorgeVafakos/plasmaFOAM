fv::options& fvOptions(fv::options::New(meshAir));

if (!fvOptions.optionList::size())
{
    Info << "No finite volume options present" << endl;
}
