C-----------------------------------------------------------------------
C  Routine:    IVOUT
C
C  Purpose:    Integer vector output routine.
C
C  Usage:      CALL IVOUT (LOUT, N, IX, IDIGIT, IFMT)
C
C  Arguments
C     N      - Length of array IX. (Input)
C     IX     - Integer array to be printed. (Input)
C     IFMT   - Format to be used in printing array IX. (Input)
C     IDIGIT - Print up to ABS(IDIGIT) decimal digits / number. (Input)
C              If IDIGIT .LT. 0, printing is done with 72 columns.
C              If IDIGIT .GT. 0, printing is done with 132 columns.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IVOUT (LOUT, N, IX, IDIGIT, IFMT)
C     ...
C     ... SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IX(*), N, IDIGIT, LOUT
      CHARACTER  IFMT*(*)
      RETURN
      END
