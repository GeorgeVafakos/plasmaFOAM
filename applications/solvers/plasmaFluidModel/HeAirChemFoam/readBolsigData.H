#include "IOdictionary.H"
#include "List.H"

Info<< "Reading bolsigProperties\n" << endl;
IOdictionary bolsigProperties
(
    IOobject
    (
        "bolsigProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

scalarList E_N_vector
(
    bolsigProperties.lookup("E_N")
);

scalarList He_vector
(
    bolsigProperties.lookup("He")
);

scalarList N2_vector
(
    bolsigProperties.lookup("N2")
);

scalarList O2_vector
(
    bolsigProperties.lookup("O2")
);

scalarList mob_vector
(
    bolsigProperties.lookup("mob")
);

// List<List<scalar>> tableData(E_N_vector.size());

// forAll(E_N_vector, i)
// {
//     tableData[i].setSize(He_vector.size());
//     forAll(He_vector, j)
//     {
//         tableData[i][j] = mob_vector[i*E_N_vector.size() + j];
//     }
// }

// // Info<< "tableData = \n" << tableData << endl;
// interpolation2DTable<scalar> mobilityTable(
//     tableData,
//     bounds::normalBounding(bounds::normalBounding::CLAMP),
//     fileName("mobilityTable")              // File name identifier (can be dummy if unused)
// );


// List<Tuple2<scalar, List<Tuple2<scalar, scalar>>>> interpTable(He_vector.size());

// forAll(He_vector, i)
// {
//     List<Tuple2<scalar, scalar>> row(E_N_vector.size());
//     forAll(E_N_vector, j)
//     {
//         // Flattened index in 1D mob_vector
//         label index = i * E_N_vector.size() + j;
//         row[j] = Tuple2<scalar, scalar>(E_N_vector[j], mob_vector[index]);
//     }

//     interpTable[i] = Tuple2<scalar, List<Tuple2<scalar, scalar>>>(He_vector[i], row);
// }

// interpolation2DTable<scalar> mobilityTable(
//     interpTable,
//     bounds::normalBounding(bounds::normalBounding::CLAMP),
//     fileName("mobilityTable")              // File name identifier (can be dummy if unused)
// );

// scalar mu_interp = mobilityTable(0.9, 96.0);

// Info<< "Interpolated mobility at E/N=90.0 and He=0.4: " << mu_interp << endl;

