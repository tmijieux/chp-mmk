#ifndef CHP_PROC_H
#define CHP_PROC_H


namespace chp {

class proc {
public:
    virtual int rank() const = 0;
    virtual int size() const = 0;
    virtual ~proc() {};
};

class mpi_proc : public proc {
    int m_rank;
    int m_group_size;

public:
    inline int rank() const override { return m_rank; }
    inline int size() const override { return m_group_size; }

    mpi_proc();
    virtual ~mpi_proc();
};

class unique_proc : public proc {
public:
    inline int rank() const override { return 0; }
    inline int size() const override { return 1; }
    unique_proc();
};

};

#endif // CHP_PROC_H
