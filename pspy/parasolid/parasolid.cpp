#include "parasolid.h"
int parasolid_session_initialized = 0;

void initialize_parasolid_frustrum()
{
    PK_SESSION_frustrum_t fru;
    PK_SESSION_frustrum_o_m(fru);
    fru.fstart = FSTART;
    fru.fstop =  FSTOP;
    fru.fmallo = FMALLO;
    fru.fmfree = FMFREE;
    fru.ffoprd = FFOPRD;
    fru.ffopwr = FFOPWR;
    fru.ffclos = FFCLOS;
    fru.ffread = FFREAD;
    fru.ffwrit = FFWRIT;
    PK_SESSION_register_frustrum(&fru);
}

int initialize_parasolid_session()
{
    initialize_parasolid_frustrum();
    PK_ERROR_t err = PK_ERROR_no_errors;
    PK_SESSION_start_o_t start_options;
    PK_SESSION_start_o_m(start_options);
    err = PK_SESSION_start(&start_options);
    return err;
}

void ensure_parasolid_session()
{
    if (!parasolid_session_initialized) {
        int err = initialize_parasolid_session();
        if (err == 0) {
            parasolid_session_initialized = 1;
        }
    }
}
