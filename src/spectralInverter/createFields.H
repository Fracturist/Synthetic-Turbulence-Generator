// Boltzmann constant 1.380649 x 10^-23
Random rndGen(1380649);

Info<< nl;

Info<< "Reading/initializing mean velocity field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector(dimVelocity, Zero)
);

Info<< "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    fvc::flux(U)
);

Info<< "Reading dynamic viscosity nu" << nl;
singlePhaseTransportModel transport(U, phi);
const scalarField& nu = transport.nu();
Info<< endl;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const word tField = args.getOrDefault<word>("turb", "k");

tmp<volScalarField> tk;
tmp<volSymmTensorField> tR;

if (tField == "k")
{
    Info<< "Use turbulent kinetic energy " << tField << nl;

    tk =
        tmp<volScalarField>::New
        (
            IOobject
            (
                "k",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

    tR =
        tmp<volSymmTensorField>::New
        (
            IOobject
            (
                "R",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedSymmTensor(dimVelocity*dimVelocity, Zero)
        );

    forAll(tk(), cell)
    {
        const scalar uu = 2.0/3.0*tk()[cell];

        tR.ref()[cell].xx() = uu;
        tR.ref()[cell].yy() = uu;
        tR.ref()[cell].zz() = uu;
    }

}
else if (tField == "R")
{
    Info<< "Use Reynolds stress tensor " << tField << nl;

    tR =
        tmp<volSymmTensorField>::New
        (
            IOobject
            (
                "R",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

    tk =
        tmp<volScalarField>::New
        (
            IOobject
            (
                "k",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimVelocity*dimVelocity, Zero)
        );

    forAll(tk(), cell)
    {
        tk.ref()[cell] = 0.5*tr(tR()[cell]);
    }
}
else
{
    FatalErrorIn(args.executable())
        << "Unknown turbulent field type " << tField
        << exit(FatalError);
}
Info<< endl;

tmp<volScalarField> tepsilon;

if (args.found("readEpsilon"))
{
    Info<< "Reading turbulent dissipation field epsilon\n" << endl;

    tepsilon =
        tmp<volScalarField>::New
        (
            IOobject
            (
                "epsilon",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const word dictName = "spectralMethodDict";
IOdictionary dict
(
    IOobject
    (
        dictName,
        runTime.constant(),
        runTime,
        IOobject::MUST_READ
    )
);

// Number of modes
const label nModes = dict.getOrDefault<label>("nModes", 5000);

// Logarithmic distribution method
const bool logKappa = dict.getOrDefault<Switch>("logDistribution", true);

// Energy spectrum form
const Switch useSpectrumTable = dict.getOrDefault("useSpectrumTable", false);

Info<< "========================================" << nl
    << "Spectral method details: " << nl
    << "  " << "Number of modes is " << nModes << nl
    << "  " << "Employing " << (logKappa ? "logarithmic" : "linear") << " distribution" << nl
    << "  " << "Energy spectrum given in " << (useSpectrumTable ? "TABLE" : "MODEL") << " form" << nl
    << "========================================" << nl << endl;

scalarField l, eta;
if (args.found("readEpsilon"))
{
    l = Foam::pow(2.0/3.0*tk(), 1.5)/tepsilon();
    eta = Foam::pow(tepsilon(), -0.25)*Foam::pow(nu, 0.75);
}

tepsilon.clear();

// Wave number vs energy profile
autoPtr<Function1<scalar>> EkTab;
autoPtr<spectrumModel> EkMod;
scalar uTur0 = 0;

if (useSpectrumTable)
{
    EkTab = Function1<scalar>::New( "Ek", dict.subDict("tabulatedSpectrum") );
    uTur0 = dict.subDict("tabulatedSpectrum").get<scalar>("uTur0");
}
else
{
    if (args.found("readEpsilon"))
    {
        EkMod = autoPtr<spectrumModel>::New(dict, l, eta);
    }
    else
    {
        EkMod = autoPtr<spectrumModel>::New(tk(), nu, dict);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tk.clear();

vectorField Uc(U.size(), Zero);

const vectorField& C(mesh.C());

Info<< "Computing transformation tensor field L" << endl;

tmp<volTensorField> tL
(
    tmp<volTensorField>::New
    (
        IOobject
        (
            "L",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor(dimVelocity, Zero)
    )
);

// Cholesky decomposition
tL.ref().primitiveFieldRef() = Chol(tR());

tL.ref().correctBoundaryConditions();

Info<< endl;

tR.clear();
