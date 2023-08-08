// Extents in x, y, z directions

const vector range(mesh.bounds().max() - mesh.bounds().min());

Info<< "range = " << range << endl;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Minimum cell spacing in x, y and z directions

vector delta(Zero);

{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    vectorField cellSpan(cells.size());

    forAll(cells, i)
    {
        const cell& c = cells[i];
        boundBox cellBb(boundBox::invertedBox);
        cellBb.add(points, c.labels(faces));
        cellSpan[i] = cellBb.span();
    }

    delta = gMin(cellSpan);
}

Info<< "delta = " << delta << endl;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Minimum wave number
scalar kappaMin = mathematical::twoPi/cmptMax(range);

// Maximum wave number
scalar kappaMax = mathematical::pi/cmptMin(delta);

Info<< "Wave number min/max = " << kappaMin << ", " << kappaMax << endl;
