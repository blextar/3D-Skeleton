#ifndef LAP_HULL_H
#define	LAP_HULL_H

namespace lap_hull{
    class face{
    public:
        double normal[3];
        double area;
        face(node *n1, node *n2, node *n3);
    };

    class vertex{
    public:
        face *f[4];
        vertex(face *f1, face *f2, face *f3, face *f4);
    };
};

#endif

