# Local relaxation of the one-dimensional Fredrickson-Andersen model for infection density above one third

Yitzchak Shmalo

In the one-spin facilitated Fredrickson-Andersen model on Z with equilibrium density q of infected sites, the law of
the process seen in a fixed finite window converges to equilibrium exponentially fast for every q > 1/3, up to an error
equal to twice the probability that no infection starts within distance vt of the origin, where v = v(q) > 0. It
follows that both conjectures of Martinelli, Shapira and Toninelli hold for q > 1/3: the process converges from every
initial law with an infection almost surely, and it converges exponentially fast when the distance from the origin to
the nearest infection has an exponentially decaying tail. Previously this was known for q > 1/2.

- `main.tex`, `main.pdf`: the paper.
- `checks/drift_check.py`: an exact rational check of the drift computation for the marked infection (Lemma 3.2 and
  Remark 6.1). It applies the rules of Definition 3.1 to every local configuration, confirms that the maximum drift is
  -q(3q-1)/(1+q) for q = j/300 with 100 < j < 300, and that no choice of the weight u gives a negative drift when
  q <= 1/3. It needs Python 3 and nothing else, and runs in a few minutes.

See `CREDITS.md` for the work this paper builds on.
