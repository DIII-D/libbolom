
        real*4 Function basisel(iparm,x,tens,rknots,numknt)

				implicit none
				integer iparm,numknt

C 20060607 tbt Changed dimension of input argument from (1) to (*)
C				real*4 x,rknots(1)
				real*4 x,rknots(*)

        integer nk
				real*4 w,tens2,tens

        basisel = 0.0
              nk = ((iparm - 1) / 2) + 1
							tens2 = abs(tens)*float(numknt-1)/
     .                     (rknots(numknt)-rknots(1)) 
							if (nk .gt. 1 .and. x .le. rknots(nk) .and. 
     .             x .ge. rknots(nk-1)) then
                    w = rknots(nk) - rknots(nk-1)
                    if(mod(iparm,2) .eq. 0) then
                          basisel = (sinh(tens2*(x-rknots(nk-1)))/
     .                          sinh(tens2*w) - (x-rknots(nk-1))/w)
     .                            / (tens2*tens2)
										else
                          basisel = (x-rknots(nk-1))/w
										endif
										
               endif
							if (nk .lt. numknt .and. x .ge. rknots(nk) .and. 
     .             x .le. rknots(nk+1)) then
                    w = rknots(nk+1) - rknots(nk)
                    if(mod(iparm,2) .eq. 0) then
                          basisel = (sinh(tens2*(rknots(nk+1)-x))/
     .                         sinh(tens2*w) - (rknots(nk+1)-x)/w)
     .                            / (tens2*tens2)
										else
                          basisel = (rknots(nk+1) - x)/w
										endif
										
               endif
        return
        end

C------------------------------------------------------------------

        real*4 Function basisv(numunk,x,tens,sol,rknots,numknt)

				implicit none

				integer numknt
				real*4 x,sol(*),rknots(*),tens
				real*4 basisel
				integer numunk,i


        basisv = 0.0
        do i = 1, numunk
             basisv = basisv+sol(i)*
     *		      basisel(i,x,tens,rknots,numknt)
        enddo
        return
        end
c
c
        
        subroutine bacnst(numunk,tens,rknots,
     *			  numknt,crsp,ncrsp,eqn,z)

				implicit none
				integer numunk,numknt,ncrsp,i,j,eqn
				real*4 rknots(*),h,w
        real*4 crsp(ncrsp ,numunk),z(ncrsp),tens,tens2
				real*4 basisel

         do j = 1,ncrsp
	   do i = 1,numunk
            crsp(j,i) = 0.0
           enddo
	   z(j) = 0.0
         enddo

	 eqn = 0
            if(numknt .le. 2)return
c
c first set of constraints is that splines have equal first 
c derivative at the rknots
c
							tens2 = abs(tens)*float(numknt-1)/
     .                     (rknots(numknt)-rknots(1)) 
            do i = 2,numknt-1
               eqn = eqn + 1
               z(eqn) = 0.0
							 w = rknots(i+1) - rknots(i)
               crsp(eqn,2*(i-1) + 1) = -1.0/w
               crsp(eqn,2*(i-1) + 2) = (-tens2*
     .   							 cosh(tens2*w)/sinh(tens2*w)
     .                 		 + 1.0/w)/(tens2*tens2)
               crsp(eqn,2*(i-1) + 3) = 1.0/w
               crsp(eqn,2*(i-1) + 4) = (tens2/
     .   							 sinh(tens2*w) - 1.0/w)/(tens2*tens2)

							 w = rknots(i) - rknots(i-1)
               crsp(eqn,2*(i-1) - 1) = 1.0/w
               crsp(eqn,2*(i-1) + 0) = -(-tens2/
     .   							 sinh(tens2*w) + 1.0/w)/(tens2*tens2)
               crsp(eqn,2*(i-1) + 1) = 
     .         crsp(eqn,2*(i-1) + 1) - 1.0/w
               crsp(eqn,2*(i-1) + 2) = 
     .               crsp(eqn,2*(i-1) + 2) - (tens2*
     .   							 cosh(tens2*w)/sinh(tens2*w)
     .             		 - 1.0/w)/(tens2*tens2)
            enddo

               eqn = eqn + 1
               z(eqn) = 0.0
							 w = rknots(2) - rknots(1)
               crsp(eqn,1) = -1.0/w
               crsp(eqn,2) = (-tens2*
     .   							 cosh(tens2*w)/sinh(tens2*w)
     .                 		 + 1.0/w)/(tens2*tens2)
               crsp(eqn,3) = 1.0/w
               crsp(eqn,4) = (tens2/
     .   							 sinh(tens2*w) - 1.0/w)/(tens2*tens2)

               eqn = eqn + 1
               z(eqn) = 0.0
							 i = numknt
							 w = rknots(numknt) - rknots(numknt - 1)
               crsp(eqn,2*(i-1) + 1) = -1.0/w
               crsp(eqn,2*(i-1) + 2) = (-tens2*
     .   							 cosh(tens2*w)/sinh(tens2*w)
     .                 		 + 1.0/w)/(tens2*tens2)
               crsp(eqn,2*(i-1) + 3) = 1.0/w
               crsp(eqn,2*(i-1) + 4) = (tens2/
     .   							 sinh(tens2*w) - 1.0/w)/(tens2*tens2)

               eqn = eqn + 1
              do i = 1, numunk
                 z(eqn) = 0.0
                  crsp(eqn,i) =  basisel(i,
     .										 rknots(numknt),tens,rknots,numknt)
              enddo
         return
         end
