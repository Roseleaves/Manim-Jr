#set page(columns:3)
#set text(6pt)
#show math.equation:it => {
  if it.fields().keys().contains("label"){
    math.equation(block: true, numbering: "(1)", it)
  } else {
    it
  }
}
$ sin x < x < tan x quad(exists a>0, 0<x<a) $  <ineq>
$ sin (alpha plus.minus beta) &= 
  sin alpha cos beta plus.minus cos alpha sin beta\
  cos (alpha plus.minus beta) &= 
  cos alpha cos beta minus.plus sin alpha sin beta\
  // tan (alpha plus.minus beta) &= 
  // (tan alpha plus.minus tan beta)/
  // (1 minus.plus tan alpha tan beta) 
$ <SCTab>
$ cos 0 = 1 quad 
  sin 0 = 0 $ <special-value>
$ cos (-beta)= cos beta quad 
  sin (-beta)=-sin beta $  <odd-even>
$ cos^2 alpha + sin^2 alpha=1 $ <sum-of-square>
$ cos alpha cos beta = 1/2[cos (alpha + beta) + cos (alpha - beta)] $ <p2s>
$ sin theta - sin phi = 2 cos (theta+phi)/2 sin (theta-phi)/2 $ <s2p>
$ tan (theta + phi)/2 = 
  (sin theta + sin phi)/(cos theta + cos phi) = 
 -(cos theta - cos phi)/(sin theta - sin phi) $ <half>
$ cos 2 alpha &= 2 cos^2 alpha-1 &quad cos n alpha &= T_n (cos alpha) \
 sin 2 alpha &= 2 sin alpha cos alpha &quad sin n alpha &= dots $ <na>
$ T_0 =1 quad T_1 = id quad T_n = 2 T_(n-1) id - T_(n-2) $ <Tn>
$ tan alpha = (2t)/(1-t^2) sin alpha = (2t)/(1+t^2) cos alpha = (1-t^2)/(1+t^2) $ <allpropose>
// $ c_0 &= 1/2, quad c_n=sqrt((1+c_(n-1))/2) -> 1 $ <circle-cutter>
$ sin' &=  cos quad & arcsin' &=  (1-id^2)^(-1/2),\
  cos' &= -sin quad & arccos' &= -(1-id^2)^(-1/2),\
  tan' &=sec^2 quad & arctan' &=1/(1+id^2) $ <deravitve>
