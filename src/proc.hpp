#ifndef CHP_PROC_H
#define CHP_PROC_H


namespace chp {

class proc {
    int m_rank;
    int m_group_size;
    
public:
    inline int rank() const { return m_rank; }
    inline int size() const { return m_group_size; } 
    
    proc();
    virtual ~proc();
};

};

#endif // CHP_PROC_H
