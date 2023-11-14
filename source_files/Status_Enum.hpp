/*
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

Header file for enumeration of 'status', which refers to whether the
simulation is currently in a 'dialling in' phase, or is not dialling and is
instead waiting for equilibrium, or whether equilibrium has been reached but
the next dialling phase has not yet begun.*/

#ifndef _STATUS_ENUM_TAG_
#define _STATUS_ENUM_TAG_

enum Status_Enum {dialling, waiting_for_equil, equil_reached};

#endif
