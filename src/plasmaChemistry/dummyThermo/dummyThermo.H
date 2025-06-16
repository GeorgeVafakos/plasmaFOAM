#ifndef dummyThermo_H
#define dummyThermo_H

#include "specie.H"

namespace Foam
{

class dummyThermo : public specie
{
public:
    TypeName("dummyThermo");

    using thermoType = dummyThermo;

    dummyThermo(const word& name);
    dummyThermo(const dictionary& dict);
    dummyThermo(const specie& sp);  // ✅ for Reaction compatibility
    virtual ~dummyThermo() {}

    autoPtr<dummyThermo> clone() const;
};

}

#endif
