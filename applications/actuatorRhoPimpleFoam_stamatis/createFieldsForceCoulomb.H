//*******************************************************

IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            meshAir,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
wordList forcesNames
(
	physicalProperties.lookup("forcesNames")
);


int counter = 0;

PtrList <volVectorField> fCoulomb(forcesNames.size());

forAll(forcesNames, s)
{
	fCoulomb.set
	(
		s,
		new volVectorField
		(
			IOobject
			(
				forcesNames[s],
				runTime.constant(),
				mesh,
				IOobject::MUST_READ,
        		IOobject::AUTO_WRITE
			),
		mesh
		)
	);
}


//*******************************************************


