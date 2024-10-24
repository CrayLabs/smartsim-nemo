MODULE ssdb
   !!======================================================================
   !!                     ***  MODULE  ssdb  ***
   !! Ocean domain  :  1D configuration
   !!=====================================================================
   !! History :  5.0  !  2024-10 (A. Shao)     Original code
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ssdb_init      : initialize SmartRedis client and other variables
   !!----------------------------------------------------------------------

   USE lib_mpp

   USE smartredis_client

   IMPLICIT NONE
   PRIVATE

   PUBLIC ssdb_init

   LOGICAL, PUBLIC :: ln_ssdb = .false.
   CHARACTER(len=16), PUBLIC :: client_prefix
   TYPE(client_type), PUBLIC :: sr_client

   CONTAINS

   SUBROUTINE ssdb_init
      INTEGER :: rc

      write(client_prefix, "(I6.6)") mpprank
      rc = sr_client%initialize()
      if (sr_client%sr_error_parser(rc)) call ctl_stop('SmartSim database client could not be initialized')

   END SUBROUTINE ssdb_init

END MODULE ssdb