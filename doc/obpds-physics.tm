<TeXmacs|1.99.2>

<style|generic>

<\body>
  <section|Carrier Density>

  <subsection|Parabolic Band Approximation>

  For a parabolic conduction-band minimum (valley) with an effective
  conduction-band density of states given by

  <\equation>
    N<rsub|C>=2<around*|(|<frac|m<rsub|e><rsup|\<ast\>>*k*T|2*\<pi\>*\<hbar\><rsup|2>>|)><rsup|3/2>,
  </equation>

  where <math|m<rsub|e><rsup|\<ast\>>> is the electron density of states
  effective mass, the electron density in the conduction-band is given by

  <\equation>
    n<around*|(|\<phi\>|)>=<frac|2*N<rsub|C>|<sqrt|\<pi\>>>F<rsub|1/2><around*|(|\<phi\>|)>.
  </equation>

  The variable, <math|\<phi\>>, in this equation is the normalized Fermi
  energy, given by

  <\equation*>
    \<phi\>=<around*|(|E<rsub|F>-E<rsub|C>|)>/k*T,
  </equation*>

  and the function <math|F<rsub|1/2>> is the Fermi-Dirac integral of order
  k=1/2, which is given by

  <\equation>
    F<rsub|1/2><around*|(|\<phi\>|)>=<big|int><rsub|0><rsup|\<infty\>><frac|\<epsilon\><rsup|1/2>*|exp<around*|(|\<epsilon\>-\<phi\>|)>+1>*d\<epsilon\>.
  </equation>

  The variable, <math|\<epsilon\>>, in this equation is the normalized
  electron kinetic energy, given by

  <\equation>
    \<epsilon\>=<around*|(|E-E<rsub|C>|)>/k*T.
  </equation>

  Similarly, in a parabolic valance-band maximum with effective valance-band
  density of states, <math|N<rsub|V>>, given by

  <\equation>
    N<rsub|V>=2<around*|(|<frac|m<rsub|h><rsup|\<ast\>>*k*T|2*\<pi\>*\<hbar\><rsup|2>>|)><rsup|3/2>,
  </equation>

  where <math|m<rsub|h><rsup|\<ast\>>> is the hole density of states
  effective mass, the hole density in the valance-band is given by

  <\equation>
    p<around*|(|\<phi\>|)>=<frac|2*N<rsub|V>|<sqrt|\<pi\>>>F<rsub|1/2><around*|(|\<phi\>|)>,
  </equation>

  where the normalized Fermi energy is given by

  <\equation*>
    \<phi\>=<around*|(|E<rsub|V>-E<rsub|F>|)>/k*T,
  </equation*>

  and the normalized hole kinetic energy is given by

  <\equation*>
    \<epsilon\>=<around*|(|E<rsub|V>-E|)>/k*T.
  </equation*>

  In OBPDS, <math|F<rsub|1/2><around*|(|\<phi\>|)>> is numerically
  approximated using the method described by Fukushima in Ref.
  <cite|fukushima_precise_2015b>.

  <subsection|Non-Parabolic Kane Approximation>

  For III-V compound semiconductors with the zinc blende crystal structure,
  e.g. GaAs, AlAs, InSb, etc., band non-parabolicity of the conduction-band
  minimum at zone center, i.e. of the \<Gamma\>-valley, can have a
  significant effect on the physics of electronic and optoelectronic devices.
  This band non-parabolicity can be modeled by the Kane k.p approximation.
  Under this approximation, the band dispersion is given by

  <\equation>
    E<rsub|k><around*|(|1+\<alpha\>*E<rsub|k>|)>=<frac|\<hbar\><rsup|2>*k<rsup|2>|2*m<rsub|e><rsup|\<ast\>>>,
  </equation>

  \;

  \ where \<alpha\> is the nonparabolicity parameter given (approximately) by

  <\equation>
    \<alpha\>=<frac|1|\<epsilon\><rsub|g>><around*|(|1-<frac|m<rsub|e><rsup|*\<ast\>>|m<rsub|e>>|)><rsup|2>,
  </equation>

  where <math|\<epsilon\><rsub|g>=E<rsub|g>/k*T>. Due to the increased
  density of states, the electron density, <math|n>, is higher than for a
  parabolic band. The relationship between the Fermi energy and the electron
  density under the Kane approximation is given by

  <\equation>
    n<around*|(|\<phi\>,\<alpha\>|)>=<frac|2*N<rsub|C>|<sqrt|\<pi\>>>H<rsub|1/2><around*|(|\<phi\>,\<alpha\>|)>,
  </equation>

  where

  <\equation>
    H<rsub|1/2><around*|(|\<phi\>,\<alpha\>|)>=<big|int><rsub|0><rsup|\<infty\>><frac|\<epsilon\><rsup|1/2>*<around*|(|1+\<alpha\>*\<epsilon\>|)><rsup|1/2><around*|(|1+2*\<alpha\>*\<epsilon\>|)>|exp<around*|(|\<epsilon\>-\<phi\>|)>+1>*d\<epsilon\>.
  </equation>

  This can expanded into

  <\equation>
    H<rsub|1/2><around*|(|\<phi\>,\<alpha\>|)>=G<rsub|1/2><around*|(|\<phi\>,\<alpha\>|)>+2*\<alpha\>*G<rsub|1><around*|(|\<phi\>,\<alpha\>|)>,
  </equation>

  where <math|G<rsub|k><around*|(|\<phi\>,\<alpha\>|)>> is a generalized
  Fermi-Dirac integral common in relativistic physics, given by

  <\equation>
    G<rsub|k><around*|(|\<phi\>,\<alpha\>|)>=<big|int><rsub|0><rsup|\<infty\>><frac|\<epsilon\><rsup|k>*<around*|(|1+\<alpha\>*\<epsilon\>|)><rsup|1/2>|exp<around*|(|\<epsilon\>-\<phi\>|)>+1>*d\<epsilon\>.
  </equation>

  In OBPDS, <math|G<rsub|k><around*|(|\<phi\>,\<alpha\>|)>> is numerically
  approximated using the method described by Fukushima in Ref.
  <cite|fukushima_precise_2015c>.

  <\bibliography|bib|tm-plain|obpds-physics.bib>
    <\bib-list|2>
      <bibitem*|1><label|bib-fukushima_precise_2015b>Toshio<nbsp>Fukushima.<newblock>
      Precise and fast computation of Fermi\UDirac integral of integer and
      half integer order by piecewise minimax rational
      approximation.<newblock> <with|font-shape|italic|Applied Mathematics
      and Computation>, 259:708--729, May 2015.<newblock>

      <bibitem*|2><label|bib-fukushima_precise_2015c>Toshio<nbsp>Fukushima.<newblock>
      Precise and fast computation of generalized Fermi\UDirac integral by
      parameter polynomial approximation.<newblock>
      <with|font-shape|italic|Applied Mathematics and Computation>,
      270:802--807, Nov 2015.<newblock>
    </bib-list>
  </bibliography>
</body>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|1.1|?>>
    <associate|auto-3|<tuple|1.2|?>>
    <associate|auto-4|<tuple|12|?>>
    <associate|auto-5|<tuple|12|?>>
    <associate|bib-fukushima_precise_2014|<tuple|1|?>>
    <associate|bib-fukushima_precise_2015b|<tuple|1|?>>
    <associate|bib-fukushima_precise_2015c|<tuple|2|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      fukushima_precise_2015b

      fukushima_precise_2015c
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Carrier
      Density> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1tab>|1.1<space|2spc>Parabolic Band Approximation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1tab>|1.2<space|2spc>Non-Parabolic Kane
      Approximation <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>