#ifndef FRUSTRUM_H_INCLUDED
#define FRUSTRUM_H_INCLUDED 1

extern void FSTART(int*);
extern void FSTOP(int*);
extern void FMALLO(int*, char**, int*);
extern void FMFREE(int*, char**, int*);
extern void FFOPRD(const int*, const int*, const char*, const int*,
    const int*, int*, int*);
extern void FFOPWR(const int*, const int*, const char*, const int*,
    const char*, const int*, int*, int*);
extern void FFCLOS(const int*, const int*, const int*, int*);
extern void FFREAD(const int*, const int*, const int*, char*, int*,
    int*);
extern void FFWRIT(const int*, const int*, const int*, const char*,
    int*);

#endif // !FRUSTRUM_H_INCLUDED
