#ifndef PARASOLID_H_INCLUDED
#define PARASOLID_H_INCLUDED 1

extern "C" {
    #include "parasolid_kernel.h"
    #include "kernel_interface.h"
}
#include "frustrum.h"

int initialize_parasolid_session();
void initialize_parasolid_frustrum();
void ensure_parasolid_session();

#endif