#lang racket
(require math/matrix)
(require math/array)
;; https://zhuanlan.zhihu.com/p/721568584
(struct problem (dest relations restrictions maximize?))
(struct stmt (l sym r) #:transparent)

(define (does-specify? probl var)
  (define (has-name? p)
    (and (pair? p) (eq? (car p) var)))
  (ormap has-name? (problem-restrictions probl)))

(define doesnt-specify? (negate does-specify?))

(define (varname? n)
  (and (symbol? n) (string-prefix? (symbol->string n) "x")))

(define (var-index n)
  (if (not (varname? n))
      0
      (string->number (substring (symbol->string n) 1))))

(define (varpair? p)
  (and (pair? p) (varname? (cdr p))))

(define (stmt-equation? st)
  (eq? '= (stmt-sym st)))

(define (non-negative-restraint? x)
  (and (pair? x) (eq? '>=0 (cdr x))))

(define (expr-sort lst)
  (sort lst
        symbol<?
        #:key (λ (x)
                (if (pair? x)
                    (cdr x)
                    'zzzzzz))))
(define (restrictions-sort lst)
  (sort lst symbol<? #:key car))
(define/match (stmt-sort st)
  [((stmt l s r)) (stmt (expr-sort l) s (expr-sort r))])

(define/match (problem-sort p)
  [((problem dests relations restrictions maximize))
   (problem (expr-sort dests) (map stmt-sort relations) (restrictions-sort restrictions) maximize)])

(define/match (standard-form? p)
  [((problem dest relations restrictions _))
   (and (andmap stmt-equation? relations)
        (andmap non-negative-restraint? restrictions)
        (andmap (curry does-specify? probl) (all-vars probl)))])

(define/match (max-var-index x)
  [((problem dest relations restrictions _)) (apply max (map max-var-index (cons dest relations)))]
  [((stmt l sym r)) (apply max (map max-var-index (append l r)))]
  [((list xs ...)) (foldl (λ (y acc) (max (max-var-index y) acc)) 0 xs)]
  [((? varpair? p)) (var-index (cdr p))]
  [((? number? _)) 0]
  [(_) 0])

(define (number->variable n)
  (string->symbol (~a "x" n)))

(define (next-variable p)
  (define n (max-var-index p))
  (number->variable (add1 n)))

(define (all-vars x)
  (define (build-vars n)
    (build-list n (λ (i) (number->variable (add1 i)))))
  (build-vars (max-var-index x)))

(define (var-to-latex v)
  (define n (var-index v))
  (format (if (< n 10) "x_~a" "x_{~a}") n))
; Min $f = -x_1 + 3x_2  s.t. \left\{ \begin{aligned}
;        4x_1 + 7x_2 \geq 56 \\
;        3x_1 + -5x_2 \geq 15\\
;        x_1 , x_2 \geq 0
;    \end{aligned}\right.$
(define (expr-to-latex expr)
  (cond
    [(varpair? expr)
     (define n (var-index (cdr expr)))
     (define effi (car expr))
     (define fmtstr (if (< n 10) "~a~ax_~a" "~a~ax_{~a}"))
     (format fmtstr
             (if (negative? effi) "" "+")
             (match effi
               [-1 "-"]
               [1 ""]
               [_ effi])
             n)]
    [(number? expr) (~a expr)]
    [(list? expr)
     (define s (apply string-append (map expr-to-latex expr)))
     (if (string-prefix? s "+")
         (substring s 1)
         s)]))
(define/match (stmt-to-latex stmts)
  [((stmt l sym r))
   (~a (expr-to-latex l)
       (match sym
         ['>= "\\geq"]
         ['<= "\\leq"]
         ['= "="])
       (expr-to-latex r))])
(define (intersperse sep lst)
  (cond
    [(null? lst) '()] ; If the list is empty, return an empty list
    [(null? (cdr lst)) lst] ; If the list has only one element, return it
    [else
     (cons (car lst) ; Add the first element
           (cons sep (intersperse sep (cdr lst))))]))
(define (restr-to-latex res)
  (~a (var-to-latex (car res))
      (match (cdr res)
        ['>=0 " \\geq 0"]
        ['<=0 " \\leq 0"])))

(define (join sep lst)
  (apply string-append (intersperse sep (map ~a lst))))

(define (to-latex x)
  (match x
    [(problem dest relations restrictions maximize)
     (define str-restrs
       (cond
         [(standard-form? x) (join ", " (map (compose var-to-latex car) restrictions))]
         [else (apply string-append (intersperse "," (map restr-to-latex restrictions)))]))
     (~a (if maximize "Max \\(z=" "Min \\(z=")
         (expr-to-latex dest)
         " , s.t. \\left\\{\\begin{aligned}\n  "
         (apply string-append (intersperse " \\\\\n  " (map stmt-to-latex relations)))
         " \\\\\n"
         str-restrs
         "\n\\end{aligned}\\right.\\)")]))

(define (print-problem p [title #f])
  (match p
    [(problem dest relations restrictions _)
     (when title
       (display "problem: ")
       (displayln title))
     (display "dest: ")
     (display (if (problem-maximize? p) "maximize " "minimize "))
     (displayln dest)
     (displayln "relations:")
     (for-each displayln relations)

     (displayln "restrictions:")
     (for-each displayln restrictions)
     (displayln "")
     (displayln (to-latex p))]))

(define/match (expr-flip-signs expr)
  [((? list? l)) (map expr-flip-signs l)]
  [((? varpair? (cons effi name))) (cons (- effi) name)]
  [(_) expr])

;; step1 最大化转最小化
(define/match (to-mimimize p)
  [((problem dest relations restrictions maximize))
   (if maximize
       (problem (expr-flip-signs dest) relations restrictions #f)
       p)])

;; step2 不等式变等式

(define/match (find-inequality-stmt p)
  [((problem dest relations restrictions maximize)) (findf (negate stmt-equation?) relations)])

(define (list-exchange l x new-x)
  (define pos (index-of l x))
  (list-set l pos new-x))

(define/match (remove-ineuqality p st)
  [((problem dest relations restrictions maximize) (stmt l (or '< '<=) r))
   (define next-var (next-variable p))
   (define next-equation (stmt (cons (cons 1 next-var) l) '= r))
   (problem dest
            (list-exchange relations st next-equation)
            (cons (cons next-var '>=0) restrictions)
            maximize)]
  [((problem dest relations restrictions maximize) (stmt l (or '> '>=) r))
   (define next-var (next-variable p))
   (define next-equation (stmt (cons (cons -1 next-var) l) '= r))
   (problem dest
            (list-exchange relations st next-equation)
            (cons (cons next-var '>=0) restrictions)
            maximize)])

(define/match (to-equalities p)
  [((problem dest relations restrictions maximize))
   (let loop ([res p])
     ;(print-problem res)
     (let ([ineql (find-inequality-stmt res)])
       (if (not ineql)
           res
           (loop (remove-ineuqality res ineql)))))])

;; step3 消除负向变量
(define/match (find-negative-variable p)
  [((problem dest relations restrictions maximize))
   (define (negative-bound? x)
     (and (pair? x) (eq? (cdr x) '<=0)))
   (define the-pair (findf negative-bound? restrictions))
   (if the-pair
       (car the-pair)
       #f)])

(define/match (remove-negative-variable p the-var)
  [((problem dest relations restrictions maximize) v)
   ;(displayln (~a the-var " is negative"))
   (define (list-substitute lst)
     (define pos (index-where lst (λ (pa) (and (pair? pa) (eq? the-var (cdr pa))))))
     (if (not pos)
         lst
         (let ([old-pair (list-ref lst pos)])
           (match old-pair
             [(cons effi var) (list-set lst pos (cons (- effi) var))]))))
   (define/match (stmt-substitute st)
     [((stmt l sym r)) (stmt (list-substitute l) sym (list-substitute r))])
   (problem (list-substitute dest)
            (map stmt-substitute relations)
            (let ([pos (index-where restrictions (λ (p) (eq? (car p) the-var)))])
              (list-set restrictions pos (cons the-var '>=0)))
            maximize)])

(define/match (eliminate-negative-variables p)
  [((problem dest relations restrictions maximize))
   (let loop ([res p])
     ;(print-problem res)
     (let ([the-var (find-negative-variable res)])
       (if (not the-var)
           res
           (loop (remove-negative-variable res the-var)))))])

;; step4 消除剩余变量
(define/match (eliminate-unbound-variables p)
  [((problem dest relations restrictions maximize))
   (let loop ([res p])
     ;(print-problem res)
     (let ([the-var (find-unbound-variable res)])
       (if (not the-var)
           res
           (loop (remove-unbound-variable res the-var)))))])
(define/match (find-unbound-variable p)
  [((problem dest relations restrictions maximize))
   (findf (λ (var) (p . doesnt-specify? . var)) (all-vars p))])

(define/match (remove-unbound-variable p the-var)
  [((problem dest relations restrictions maximize) v)
   ; (displayln (~a the-var " is unbound"))
   (define comp-var (next-variable p))
   (define (list-substitute lst)
     (define pos (index-where lst (λ (pa) (and (pair? pa) (eq? the-var (cdr pa))))))
     (if (not pos)
         lst
         (match (list-ref lst pos)
           [(cons effi var) (cons (cons (- effi) comp-var) lst)])))
   (define/match (stmt-substitute st)
     [((stmt l sym r)) (stmt (list-substitute l) sym (list-substitute r))])

   (problem (list-substitute dest)
            (map stmt-substitute relations)
            (list* (cons the-var '>=0) (cons comp-var '>=0) restrictions)
            maximize)])

(define/match (problem-dest-gradient p)
  [((problem dest relations restrictions maximize))
   (define n-vars (max-var-index p))
   (define (gradient-of-var v-num)
     (define var (number->variable v-num))
     (match (findf (λ (pa) (eq? var (cdr pa))) dest)
       [(cons effi v) effi]
       [#f 0]))
   (vector->array (build-vector n-vars (λ (i) (gradient-of-var (add1 i)))))])

(struct solution (varnum bases value))

(define/match (solution-assignments sol)
  [((solution varnum bases value))
   (for/list ([var-n (range 1 (add1 varnum))]
              [v value])
     (cons (number->variable varnum) v))])

(define (apply-solution sol probl)
  (apply-assignments (solution-assignments sol) (problem-dest probl)))

(define (apply-assignments assi x)
  (match x
    [(stmt l sym r)
     (define vl (apply-assignments l))
     (define vr (apply-assignments r))
     (match sym
       ['>= (>= l r)]
       ['<= (<= l r)]
       ['> (> l r)]
       ['< (< l r)]
       ['= (= l r)])]
    [(? list? l) (foldl + 0 (map (curryr apply-assignments assi) l))]
    [(? varpair? (cons effi var)) (or (assoc var assi) 0)]
    [(? number? n) n]))
(define/match (solution-lookup sol var)
  [((solution varnum bases value) v)
   (define n (varnum v))
   (array-ref value (vector n))])

(define/match (solution-unbases sol)
  [((solution varnum bases value)) (remove* bases (range 1 (add1 varnum)))])
(define/match (solution-exchanges sol)
  [((solution varnum bases value))
   (for*/list ([from (in-list bases)]
               [to (in-list (solution-unbases sol))])
     (cons from to))])

(define/match (satisfy? sol probl)
  [((solution varnum bases value) (stmt _ _ _)) (apply-assignments (solution-assignments sol) stmt)]
  [((solution varnum bases value) (problem dest relations restrictions maximize))
   (define/match (restr-satisfy? r)
     [((cons var '>=0)) (not (negative? (solution-lookup sol var)))]
     [((cons var '<=0)) (not (positive? (solution-lookup sol var)))])
   (and (andmap (curry satisfy? sol) relations) (andmap (restr-satisfy? restrictions)))])
;;;;;;;;;;;;;;;;;;;;;;;;;
(define s (stmt '((1 . x1) (2 . x2)) '< '(1)))

(println s)
(println (max-var-index s))

#;(define probl
    (problem '((3 . x1) (-2 . x2))
             (list (stmt '((1 . x1) (1 . x2)) '<= '(1)) (stmt '((1 . x1) (2 . x2)) '>= '(4)))
             '((x1 . >=0))
             'maximize))
(define probl
  (problem '((2 . x1) (-1 . x2) (2 . x3))
           (list (stmt '((-1 . x1) (1 . x2) (1 . x3)) '= '(4))
                 (stmt '((-1 . x1) (1 . x2) (-1 . x3)) '<= '(6)))
           '((x1 . >=0) (x2 . <=0))
           'maximize))
;(print-problem probl)
; (println (probl . doesnt-specify? . 'x1))
; (println (probl . doesnt-specify? . 'x3))
; (println (all-vars probl))
; (displayln (non-negative-restraint? '(x1 . >=0)))
; (displayln (standard-form? probl))

(define (to-standard p [title #f])
  (when title
    (displayln title))

  (define p1 (to-mimimize p))
  (when (not (eq? p p1))
    (displayln "将目标转为最小化: \n")
    (displayln (to-latex (problem-sort p1)))
    (displayln ""))
  (define p2 (to-equalities p1))
  (when (not (eq? p2 p1))
    (displayln "不等式化等式: \n")
    (displayln (to-latex (problem-sort p2)))
    (displayln ""))
  (define p3 (eliminate-unbound-variables p2))
  (when (not (eq? p3 p2))
    (displayln "去除无界变量: \n")
    (displayln (to-latex (problem-sort p3)))
    (displayln ""))
  (define p4 (eliminate-negative-variables (problem-sort p3)))
  (when (not (eq? p4 p3))
    (displayln "去除负约束变量: \n")
    (displayln (to-latex p4))
    (displayln ""))
  p4)

(define probl4-1
  (problem '((2 . x1) (-1 . x2) (2 . x3))
           (list (stmt '((-1 . x1) (1 . x2) (1 . x3)) '= '(4))
                 (stmt '((-1 . x1) (1 . x2) (-1 . x3)) '<= '(6)))
           '((x1 . <=0) (x2 . >=0))
           #f))
(to-standard probl4-1)

(define probl4-2
  (problem '((2 . x1) (1 . x2) (3 . x3) (1 . x4))
           (list (stmt '((1 . x1) (1 . x2) (1 . x3) (1 . x4)) '<= '(7))
                 (stmt '((2 . x1) (-3 . x2) (5 . x3)) '= '(-8))
                 (stmt '((1 . x1) (-2 . x3) (2 . x4)) '>= '(-1)))
           '((x1 . >=0) (x2 . <=0) (x3 . >=0))
           #f))

(to-standard probl4-2 "4-2:")

(define probl9-1
  (problem '((6 . x1) (14 . x2) (13 . x3))
           (list (stmt '((1 . x1) (4 . x2) (2 . x3)) '<= '(48))
                 (stmt '((1 . x1) (2 . x2) (4 . x3)) '<= '(60)))
           '((x1 . >=0) (x2 . >=0) (x3 . >=0))
           'maximize))

(define probl9-1-st (to-standard probl9-1 "9-1:"))
(displayln (problem-dest-gradient probl9-1-st))

(define sol (solution 4 '(1 2) (array #[0 0 1])))
(displayln (solution-unbases sol))
(displayln (solution-exchanges sol))

(define probl9-2
  (problem '((3 . x1) (-2 . x2) (-4 . x3))
           (list (stmt '((4 . x1) (5 . x2) (-2 . x3)) '<= '(22))
                 (stmt '((1 . x1) (-2 . x2) (1 . x3)) '<= '(30)))
           '((x1 . >=0) (x2 . >=0) (x3 . >=0))
           'maximize))

(define probl9-2-st (to-standard probl9-2 "9-2:"))

(define probl9-3
  (problem '((1 . x1) (1 . x2) (1 . x3))
           (list (stmt '((1 . x1) (-2 . x3)) '<= '(5))
                 (stmt '((2 . x1) (-3 . x2) (1 . x3)) '<= '(3))
                 (stmt '((2 . x1) (-5 . x2) (6 . x3)) '<= '(5)))
           '((x1 . >=0) (x2 . >=0) (x3 . >=0))
           'maximize))

(define probl9-3-st (to-standard probl9-3 "9-3:"))

(define probl10
  (problem '((4 . x1) (2 . x2) (8 . x3))
           (list (stmt '((2 . x1) (-1 . x2) (3 . x3)) '<= '(30))
                 (stmt '((1 . x1) (2 . x2) (4 . x3)) '= '(40)))
           '((x1 . >=0) (x2 . >=0) (x3 . >=0))
           'maximize))

(define probl10-st (to-standard probl10 "probl-10:"))
