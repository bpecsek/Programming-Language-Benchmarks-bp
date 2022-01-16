(declaim (optimize speed (safety 0) (debug 0)))
(sb-int:set-floating-point-modes :traps (list :divide-by-zero))

(eval-when (:compile-toplevel :load-toplevel :execute)  
  (require "sb-md5")
  (ql:quickload '(:sb-simd :serapeum) :silent t)
  (defpackage :v8
    (:use #:common-lisp)
    (:import-from #:sb-simd-avx2 #:f64 #:f64.4+ #:f64.4- #:f64.4* #:f64.4-aref
                  #:f64.4-values #:f64 #:f64.4-permute128 #:f64.4-permute)
    (:import-from #:serapeum #:->)
    (:export #:make-f64v8 #:f64v8+ #:f64v8- #:f64v8* #:f64v8))
  
  (defpackage :mandelbrot
    (:nicknames :mb)
    (:use :cl :sb-md5 :v8 :sb-simd-avx)
    (:import-from #:cl-user #:define-alien-routine
                  #:long #:int)))

(in-package #:v8)
(deftype f64v8 () '(simple-array f64 (8)))

(declaim (ftype (function (&key (:initial-element f64)
                                (:initial-contents list)) f64v8) make-f64v8)
         (inline make-f64v8))
(defun make-f64v8 (&key initial-element initial-contents)
  (cond (initial-element  (make-array 8 :element-type 'f64
                                        :initial-element initial-element))
        (initial-contents (make-array 8 :element-type 'f64
                                        :initial-contents initial-contents))
        (t (make-array 8 :element-type 'f64))))

(macrolet ((define (name inst)
              `(progn (declaim (ftype (function (f64v8 f64v8 f64v8) f64v8) ,name)
                               (inline ,name))
                      (defun ,name (u v r)
                        (declare (optimize speed (safety 0) (debug 0))
                                 (type f64v8 u v r))
                        (setf (f64.4-aref r 0) (,inst (f64.4-aref u 0) (f64.4-aref v 0))
                              (f64.4-aref r 4) (,inst (f64.4-aref u 4) (f64.4-aref v 4)))
                        r))))
  (define f64v8+ f64.4+)
  (define f64v8- f64.4-)
  (define f64v8* f64.4*))

(in-package :mandelbrot)

(defconstant +MAX-ITER+ 50)
(defconstant +VLEN+ 8)

(eval-when (:load-toplevel :compile-toplevel :execute)
  (defun get-thread-count ()
    (progn (define-alien-routine sysconf long (name int))
           (sysconf 84)))
  (defconstant +workers+ (get-thread-count)))

(declaim (ftype (function (f64v8 f64v8) u8) mbrot8)
         (inline mbrot8))
(defun mbrot8 (cr ci)
  (declare (type f64v8 cr ci))
  (let ((zr (make-f64v8))
        (zi (make-f64v8))
        (tr (make-f64v8))
        (ti (make-f64v8))
        (absz (make-f64v8)))
    (declare (type f64v8 zr zi tr ti absz))
    (loop repeat (/ +MAX-ITER+ 5)
          do (loop repeat 5
                   with tmp of-type f64v8 = (make-f64v8)
                   do (f64v8+ zr zr tmp)
                      (f64v8* tmp zi tmp)
                      (f64v8+ tmp ci zi)
                      (f64v8- tr ti tmp)
                      (f64v8+ tmp cr zr)
                      (f64v8* zr zr tr)
                      (f64v8* zi zi ti))
             (f64v8+ tr ti absz)
             (let ((terminate t))
               (loop for i of-type u8 below +VLEN+
                     when (<= (aref absz i) 4d0)
                       do (setf terminate nil)
                          (loop-finish))
               (if terminate (return-from mbrot8 0))))
    (loop with accu of-type u8 = (u8 0)
          for i of-type u8 below +VLEN+
          do (setf accu (logior accu (if (<= (aref absz i) 4d0) (ash #x80 (- i)) 0)))
          finally (return accu))))

(declaim (ftype (function ((simple-array f64 (*)) u32) f64v8) array-slice)
         (inline array-slice))
(defun array-slice (array row)
  (loop with new-array of-type f64v8 = (make-f64v8)
        for i of-type u32 from (* row 8) below (+ (* row 8) 8)
        for j of-type u32 below 8
        do (setf (f64-aref new-array j) (f64-aref array i))
        finally (return (the f64v8 new-array))))

(defun main (n-input)
  (let* ((size (* (floor (+ n-input 7) +VLEN+) +VLEN+)) ;Round sizes to multiple of 8
         (chunk-size (floor size +VLEN+))
         (inv (/ 2d0 (f64 size)))
         (xloc (make-array size :element-type 'f64))
         (bitmap (make-array (* size chunk-size) :element-type 'u8
                                                 :initial-element (u8 0))))
    (declare (type u32 size chunk-size)
             (type f64 inv)
             (type f64vec xloc) 
             (type u8vec bitmap))
    (loop for i below size
          do (setf (aref xloc i) (- (* (f64 i) inv) 1.5d0)))
    (mapc #'sb-thread:join-thread
          (loop for i from 0 below +workers+
                collecting (sb-thread:make-thread
                            (let ((thread-num i))
                              (lambda ()
                                (loop for chunk-id of-type u32
                                      from thread-num below size by +workers+ 
                                      do (loop with ci of-type f64v8 =
                                                 (make-f64v8 :initial-element (- (* (f64 chunk-id) inv) 1d0))
                                               for i of-type u32 below chunk-size
                                               for r of-type u8 = (mbrot8 (array-slice xloc i) ci)
                                               when (> r 0)
                                                 do (setf (aref bitmap (+ (* chunk-id chunk-size) i)) r))))))))
    (format t "P4~%~d ~d~%" size size)
    (format t "~(~{~2,'0X~}~)~%" (coerce (sb-md5:md5sum-sequence bitmap) 'list))))

(in-package :cl-user)
(defun main (&optional n-supplied)
  (let* ((args sb-ext:*posix-argv*)
         (n-input (or n-supplied (parse-integer (or (car (last args)) "8000")))))
    (mandelbrot::main n-input)))

#+(or)
(let ((a (v8:make-f64v8 :initial-contents '(1d0 2d0 3d0 4d0 5d0 6d0 7d0 8d0)))
      (b (v8:make-f64v8 :initial-contents '(2d0 3d0 4d0 5d0 6d0 7d0 8d0 9d0)))
      (r (v8:make-f64v8)))
  (declare (optimize speed (safety 0)))
  (benchmark:with-timing (1) (loop repeat 1000000000 do (v8:f64v8+ a b r))))

#+(or)
(loop with inv = (/ 2d0 (f64 16))
      for i below 16
      for (f r) = (multiple-value-list (floor i 8))
      do (setf (aref (aref xloc f) r) (- (* (f64 i) inv) 1.5d0)))
