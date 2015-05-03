#!/usr/bin/csi
;;
;; Chicken Scheme bindings for the LAPACK routines in the ATLAS
;; library.
;;
;; Copyright 2015 Jeremy Steward
;;
;; This program is free software: you can redistribute it and/or
;; modify it under the terms of the GNU General Public License as
;; published by the Free Software Foundation, either version 3 of the
;; License, or (at your option) any later version.
;;
;; This program is distributed in the hope that it will be useful, but
;; WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;; General Public License for more details.
;;
;; A full copy of the GPL license can be found at
;; <http://www.gnu.org/licenses/>.
;;

(use srfi-4
     lapack-extras
     test)

(test-group "Test (s|d|c|z)geev_ procedures"
  (test "sgeev_ eigenvalues / loadings of identity are identity"
        (values '#f32(1 1 1)
                '#f32(0 0 0)
                '#f32(1 0 0
                      0 1 0
                      0 0 1))
        (let* ([jobvl "N"]
               [jobvr "V"]
               [N 3]
               [A (f32vector 1 0 0 0 1 0 0 0 1)]
               [WR (make-f32vector N 0)]
               [WI (make-f32vector N 0)]
               [LDVL N]
               [VL (make-f32vector (* LDVL N) 0)]
               [LDVR N]
               [VR (make-f32vector (* LDVR N) 0)]
               [LWORK (* 4 N)]
               [WORK (make-f32vector LWORK 0)]
               [INFO (s32vector 0)])
          (let-values ([(a wr wi vl vr work)
                       (sgeev jobvl
                              jobvr
                              N
                              A
                              WR
                              WI
                              VL
                              VR
                              WORK
                              LWORK
                              INFO)])
            (values wr wi vr))))

  (test "dgeev_ eigenvalues / loadings of identity are identity"
        (values '#f64(1 1 1)
                '#f64(0 0 0)
                '#f64(1 0 0
                      0 1 0
                      0 0 1))
        (let* ([jobvl "N"]
               [jobvr "V"]
               [N 3]
               [A (f64vector 1 0 0 0 1 0 0 0 1)]
               [WR (make-f64vector N 0)]
               [WI (make-f64vector N 0)]
               [LDVL N]
               [VL (make-f64vector (* LDVL N) 0)]
               [LDVR N]
               [VR (make-f64vector (* LDVR N) 0)]
               [LWORK (* 4 N)]
               [WORK (make-f64vector LWORK 0)]
               [INFO (s32vector 0)])
          (let-values ([(a wr wi vl vr work)
                        (dgeev jobvl
                               jobvr
                               N
                               A
                               WR
                               WI
                               VL
                               VR
                               WORK
                               LWORK
                               INFO)])
            (values wr wi vr))))
  )
