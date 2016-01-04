// isTracking.hh
// central place to turn on/off TRACKING switch
#define TRACKING

#define BUFFER_COLUMN_LEN 22     // room for start point and end point of each step
                                 // plus time/arclength entries for each.
#define ENDPOINT_BASE_INDEX 11
#define POSITION_SLOT 2
#define MOMENTUM_SLOT 5
#define RHS_SLOT 8

// Uncomment to test number used function calls supposing interpolation methods in
// DistChord() were available:
//#define NO_COUNT_FUNCTION_CALLS_FROM_DISTCHORD
