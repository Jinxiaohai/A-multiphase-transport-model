      FUNCTION asinh(x)
      SAVE
      if(x.gt.0) then
         ASINH=alog(x+sqrt(x**2+1.))
      else
c     a la suggestion de YP Liu:
         ASINH=-alog(-x+sqrt(x**2+1.))
      endif
      return
      end
