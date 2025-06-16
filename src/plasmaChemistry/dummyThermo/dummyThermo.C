#include "dummyThermo.H"

namespace Foam
{

defineTypeNameAndDebug(dummyThermo, 0);

dummyThermo::dummyThermo(const word& name)
:
    specie(name, 1.0, 1.0)
{}

dummyThermo::dummyThermo(const dictionary& dict)
:
    specie(dict)
{}

dummyThermo::dummyThermo(const specie& sp)
:
    specie(sp)
{}

autoPtr<dummyThermo> dummyThermo::clone() const
{
    return autoPtr<dummyThermo>(new dummyThermo(*this));
}

}
