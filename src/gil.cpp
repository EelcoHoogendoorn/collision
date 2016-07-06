
#pragma once

#include <Python.h>


class releaseGIL{
public:
    inline releaseGIL(){
        save_state = PyEval_SaveThread();
    }

    inline ~releaseGIL(){
        PyEval_RestoreThread(save_state);
    }
private:
    PyThreadState *save_state;
};


class AcquireGIL
{
public:
    inline AcquireGIL(){
        state = PyGILState_Ensure();
    }

    inline ~AcquireGIL(){
        PyGILState_Release(state);
    }
private:
    PyGILState_STATE state;
};
