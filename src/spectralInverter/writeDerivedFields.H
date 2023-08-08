{
    Info<< "Generating resolved kinetic energy field" << endl;
    volScalarField resolvedK("resolvedK", 0.5*magSqr(U));
    resolvedK.write();
    Info<< "min/max/average resolvedK = "
        << gMin(resolvedK) << ", " << gMax(resolvedK) << ", " << gAverage(resolvedK)
        << endl;
}

{
    Info<< "Generating div(U) field" << endl;
    volScalarField divU(fvc::div(U));
    divU.write();
    Info<< "min/max/average div(U) = "
        << gMin(divU) << ", " << gMax(divU) << ", " << gAverage(divU)
        << endl;
}
