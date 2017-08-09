#ifndef SEMI_METHOD_H
#define SEMI_METHOD_H

namespace Semi {

class semi_method {
public:
    virtual ~semi_method() { }
    virtual double get_error() = 0;
    virtual double get_energy() = 0;
    virtual void take_step() = 0;
    virtual void final_print() = 0;
};

} // namespace Semi

#endif // SEMI_METHOD_H
