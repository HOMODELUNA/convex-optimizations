#lang racket
(require math/matrix)
(require math)

(define (intersperse sep lst)
  (cond
    [(null? lst) '()] ; If the list is empty, return an empty list
    [(null? (cdr lst)) lst] ; If the list has only one element, return it
    [else
     (cons (car lst) ; Add the first element
           (cons sep (intersperse sep (cdr lst))))]))
(define (join sep lst)
  (apply string-append (intersperse sep (map ~a lst))))

(define (matrix-add-row m ir-from ir-to [n 1])
  (define r-from (matrix-row m ir-from))
  (define r-to (matrix-row m ir-to))
  (define new-r-to (matrix+ r-to (matrix-scale r-from n)))
  (matrix-set-row m ir-to new-r-to))

(define (matrix-scale-row m ir [k 1])
  (define r (matrix-row m ir))
  (define new-r (matrix-scale r k))
  (matrix-set-row m ir new-r))

(define (matrix-unify-row m ir)
  (define r (matrix-row m ir))
  (define k
    (for/first ([x (in-array r)]
                #:when (not (zero? x)))
      x))
  (if (or (not k) (= 1 k))
      m
      (matrix-scale-row m ir (/ 1 k))))

(define (matrix-to-latex m)
  (define (print-num x)
    (cond
      [(zero? x) " "]
      [(and (not (integer? x)) (rational? x))
       (if (negative? x)
           (format "-\\frac{~v}{~v}" (- (numerator x)) (denominator x))
           (format "\\frac{~v}{~v}" (numerator x) (denominator x)))]
      [else x]))
  (define (print-row r)
    (join " & " (map print-num r)))
  (~a "\\begin{bmatrix}\n" (join " \\\\\n" (map print-row (matrix->list* m))) "\n\\end{bmatrix}"))

(define (d-latex m)
  (displayln (matrix-to-latex m)))

(define (last-col-all->=0? m)
  (define-values (nrowls ncols) (matrix-shape m))
  (define last-col (matrix-col m (sub1 ncols)))
  (for/and ([x (in-array last-col)])
    (not (negative? x))))

(define (pivots->orders pivots n)
  (define rests (remove* pivots (range n)))
  (append pivots rests))

(define (matrix-reorder-cols m new-order)
  (define (assign-one-col m1 col-from col-to)
    (matrix-set-col m1 col-to (matrix-col m col-from)))
  (for/fold ([acc m])
            ([col-from (in-range (length new-order))]
             [col-to (in-list new-order)])
    (assign-one-col acc col-from col-to)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define m (matrix [[1 2 3 -1 1 4] [-1 1 -3 2 -1 -5]]))
(define m1 (matrix-row-echelon m))
(displayln m)
(displayln m1)
(define m2 (matrix-add-row m1 1 0 -2/3))
(displayln m2)
(define m3 (matrix-scale-row m2 1 1/3))
(displayln m3)
(define m4 (matrix-unify-row m2 1))
(displayln m4)
(displayln (matrix-to-latex m4))

(let ([a 2])
  (define m (matrix [[1 -1 2 3 1 4] [-1 2 1 -3 -1 5]]))
  (displayln m)
  (define m1 (matrix-row-echelon m))
  (displayln m1)
  (define m2 (matrix-add-row m1 1 0))
  (displayln m2)
  (d-latex m2))

(let ([a 2])
  (define m (matrix [[2 3 1 -1 1 4] [1 -3 -1 2 -1 5]]))
  (displayln m)
  (define m1 (matrix-row-echelon m #t #t 'first))
  (displayln m1)
  (d-latex m1))

(let ([a 2])
  (define m (matrix [[5 -2 4 1 0 22] [-2 1 1 0 1 30]]))
  (displayln m)
  (define m1 (matrix-row-echelon m #t #t 'first))
  (displayln m1)
  (d-latex m1))


(displayln "9-c\n")
(let ([a 2])
  (define m (matrix [[1 0 -2 1 0 0 5] [2 -3 1 0 1 0 3] [2 -5 6 0 0 1 5]]))
  (displayln m)
  (define-values (nrows ncols) (matrix-shape m))
  (for ([pivots (in-combinations (range (sub1 ncols)) nrows)])
    (define new-order (pivots->orders pivots (sub1 ncols)))
    (define m1 (matrix-reorder-cols m new-order))
    (define m2 (matrix-row-echelon m1 #t #t 'first))
    (printf "  ~a satisfy? ~a~%" pivots (last-col-all->=0? m2))
    (printf "    m1=~a~%" m1)
    (printf "    m2=~a~%" m2)
    (when (last-col-all->=0? m2)
      (d-latex m2))))


(displayln "10-1\n")
(let ([a 2])
  (define m (matrix [[2 -1 3 1 0 30] [1 2 4 0 1 40]]))
  (displayln m)
  (define-values (nrows ncols) (matrix-shape m))
  (for ([pivots (in-combinations (range (sub1 ncols)) nrows)])
    (define new-order (pivots->orders pivots (sub1 ncols)))
    (define m1 (matrix-reorder-cols m new-order))
    (define m2 (matrix-row-echelon m1 #t #t 'first))
    (printf "  ~a satisfy? ~a~%" pivots (last-col-all->=0? m2))
    (printf "    m1=~a~%" m1)
    (printf "    m2=~a~%" m2)
    (when (last-col-all->=0? m2)
      (d-latex m2))))