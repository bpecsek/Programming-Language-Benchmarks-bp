;;; Based on 2.zig
;;; Converted to Common Lisp by Bela Pecsek - 2021-12-04
(declaim (optimize (speed 3) (safety 0) (debug 0)))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (ql:quickload :sb-simd :silent t)
  (use-package :sb-simd-avx2))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (declaim (inline pos vel mass (setf pos) (setf vel) (setf mass)))
  (defstruct (body (:conc-name nil)
                   (:constructor make-body (pos vel mass)))
   (pos  (f64.4 0d0) :type f64.4)
   (vel  (f64.4 0d0) :type f64.4)
   (mass (f64.4 0d0) :type f64.4))

  (declaim (ftype (function (f64 f64 f64 f64 f64 f64 f64) body) create-body)
           (inline create-body))
  (defun create-body (x y z vx vy vz mass)
    (make-body (make-f64.4 x y z 0d0)
               (make-f64.4 vx vy vz 0d0)
               (f64.4 mass)))

  (declaim (ftype (function (body) f64) mass1)
           (inline mass1))
  (defun mass1 (body)
    (multiple-value-bind (a) (f64.4-values (mass body)) a))
  

  (defconstant +DAYS-PER-YEAR+ 365.24d0)
  (defconstant +SOLAR-MASS+ (* 4d0 pi pi))

  (defparameter *jupiter*
    (create-body 4.84143144246472090d0
                 -1.16032004402742839d0
                 -1.03622044471123109d-1
                 (* 1.66007664274403694d-3 +days-per-year+)
                 (* 7.69901118419740425d-3 +days-per-year+)
                 (* -6.90460016972063023d-5  +days-per-year+)
                 (* 9.54791938424326609d-4 +solar-mass+)))

  (defparameter *saturn*
    (create-body 8.34336671824457987d0
                 4.12479856412430479d0
                 -4.03523417114321381d-1
                 (* -2.76742510726862411d-3 +days-per-year+)
                 (* 4.99852801234917238d-3 +days-per-year+)
                 (* 2.30417297573763929d-5 +days-per-year+)
                 (* 2.85885980666130812d-4 +solar-mass+)))

  (defparameter *uranus*
    (create-body 1.28943695621391310d1
                 -1.51111514016986312d1
                 -2.23307578892655734d-1
                 (* 2.96460137564761618d-03 +days-per-year+)
                 (* 2.37847173959480950d-03 +days-per-year+)
                 (* -2.96589568540237556d-05 +days-per-year+)
                 (* 4.36624404335156298d-05 +solar-mass+)))

  (defparameter *neptune*
    (create-body 1.53796971148509165d+01
                 -2.59193146099879641d+01
                 1.79258772950371181d-01
                 (* 2.68067772490389322d-03 +days-per-year+)
                 (* 1.62824170038242295d-03 +days-per-year+)
                 (* -9.51592254519715870d-05 +days-per-year+)
                 (* 5.15138902046611451d-05 +solar-mass+)))

  (defparameter *sun* (create-body 0.0d0 0.0d0 0.0d0 0.0d0 0.0d0 0d0
			           +solar-mass+))
  (defparameter *system* (list *sun* *jupiter* *saturn* *uranus* *neptune*)))

;; Some helper functions
(declaim (ftype (function (f64.4 f64.4) f64) dot)
         (inline dot scale length-sq length_))
(defun dot (a b)
  (f64.4-hsum (f64.4* a b)))

(declaim (ftype (function (f64.4) f64) length-sq))
(defun length-sq (a)
  (dot a a))

(declaim (ftype (function (f64.4) f64) length_))
(defun length_ (a)
  (sqrt (length-sq a)))

;; Calculate the momentum of each body and conserve momentum of the system by
;; adding to the Sun's velocity the appropriate opposite velocity needed in
;; order to offset that body's momentum.
(declaim (ftype (function (list) null) offset-momentum))
(defun offset-momentum (system)
  (loop with pos = (f64.4 0)
        with sun = (car system)
        for bi in system
        do (f64.4-incf pos (f64.4* (vel bi) (mass bi)))
           (setf (vel sun) (f64.4* pos (/ (- +SOLAR-MASS+))))))

;; Advances with timestem dt = 1.0d0
;; Advance all the bodies in the system by one timestep. Calculate the
;; interactions between all the bodies, update each body's velocity based on
;; those interactions, and update each body's position by the distance it
;; travels in a timestep of 1.0d0 at it's updated velocity.
(declaim (ftype (function (list u32) null) advance))
(defun advance (system n)
  (loop repeat n do
    (loop for (bi . rest) on system do
      (dolist (bj rest)
        (let* ((pd  (f64.4- (pos bi) (pos bj)))
               (dsq (f64.4  (length-sq pd)))
               (dst (f64.4-sqrt dsq))
               (mag (f64.4/ (f64.4* dsq dst)))
               (pd-mag (f64.4* pd mag)))
          (f64.4-decf (vel bi) (f64.4* pd-mag (mass bj)))
          (f64.4-incf (vel bj) (f64.4* pd-mag (mass bi))))))
    (loop for b in system do
      (f64.4-incf (pos b) (vel b)))))

;; Output the total energy of the system.
(declaim (ftype (function (list) null) energy))
(defun energy (system)
    (let ((e 0d0))
      (declare (type f64 e))
      (loop for (bi . rest) on system do
	;; Add the kinetic energy for each body.
        (incf e (f64* 0.5d0 (mass1 bi) (length-sq (vel bi))))
        (dolist (bj rest)
	  ;; Add the potential energy between this body and every other bodies
          (decf e (f64/ (f64* (mass1 bi) (mass1 bj))
                        (length_ (f64.4- (pos bi) (pos bj)))))))
      (format t "~,9f~%" e))) ;; Output the total energy of the system.


(defconstant +DT+ 0.01d0)
(defconstant +RECIP-DT+ (/ +dt+))
;; Rescale certain properties of bodies. That allows doing
;; consequential advances as if dt were equal to 1.0d0.
;; When all advances done, rescale bodies back to obtain correct energy.
(declaim (ftype (function (list f64) null) scale-bodies))
(defun scale-bodies (system scale)
  (dolist (b system)
    (setf (mass b) (f64.4* (mass b) (f64* scale scale)))
    (setf (vel b)  (f64.4* (vel b) scale))))

(defun nbody (n)
  (let ((system *system*))
    (offset-momentum system)         
    (energy system)                     ;; Output initial energy of the system
    (scale-bodies system +DT+)          ;; Scale bodies to use time step of 1.0d0
    (advance system n)                  ;; Advance system n times
    (scale-bodies system +RECIP-DT+)    ;; Rescale bodies back to original
    (energy system)))                   ;; Output final energy of the system

(defun main (&optional n-supplied)
  (let ((n (or n-supplied (parse-integer (or (car (last #+sbcl sb-ext:*posix-argv*))
                                             "10000")))))
    (nbody n)))
