#ifndef GUROBIMODEL_H
#define GUROBIMODEL_H

class GurobiModel
{
public:
    GurobiModel();

    double objective() const;

    size_t numConstraints() const;
};

#endif  // GUROBIMODEL_H
