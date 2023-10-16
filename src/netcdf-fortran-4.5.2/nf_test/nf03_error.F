C     Copyright 1996, UCAR/Unidata
C     See netcdf/COPYRIGHT file for copying and redistribution conditions.

C     Steve Emmerson, Ed Hartnett

C
C     Use for logging error messages
C
      subroutine error(msg)
      USE tests
      implicit none
      character*(*)   msg

      nfails = nfails + 1
      if (nfails .le. max_nmpt) print *, msg
      end


C
C     Use for logging error conditions
C
      subroutine errori(msg, i)
      USE tests
      implicit none
      character*(*)   msg
      integer         i

      nfails = nfails + 1
      if (nfails .le. max_nmpt) print *, msg, i
      end


C
C     Use for logging error conditions
C
      subroutine errord(msg, d)
      USE tests
      implicit none
      character*(*)   msg
      doubleprecision d

      nfails = nfails + 1
      if (nfails .le. max_nmpt) print *, msg, d
      end


C
C     Use for logging error conditions
C
      subroutine errorc(msg, string)
      USE tests
      implicit none
      character*(*)   msg
      character*(*)   string

      nfails = nfails + 1
      if (nfails .le. max_nmpt) print *, msg,
     +     string(1:len_trim(string))
      end


C
C     Use for logging error conditions
C
      subroutine errore(msg, err)
      USE tests
      implicit none
      character*(*)   msg
      integer         err

      call errorc(msg, nf_strerror(err))
      end
