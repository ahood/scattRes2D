all: manual.pdf apply_Tstar.pdf

SECS = intro.tex mesh.tex potential.tex scattering.tex \
	scattResComp2d.tex dirBC.tex pmlBC.tex DtNBC.tex \
	scattering_sparse.tex scattResComp2d_sparse.tex \
	dirBC_sparse.tex pmlBC_sparse.tex DtNBC_sparse.tex \
	resonances.tex shapes.tex rat_approx.tex \
	ratApprox.tex \
	hankel_zeros.tex ratApproxDtNBC.tex ratApproxDtNBC_sparse.tex \
	halfcheb.tex laplacian.tex piecewise.tex dtn_map.tex \
	construct_dtnC.tex preconditioners.tex construct_RHS.tex \
	pml.tex axisymm.tex cheb_basis.tex
EPS = $(wildcard figures/*.eps)
FIGS = $(EPS:.eps=.pdf)

manual.pdf: manual.tex $(SECS) $(FIGS) 
	pdflatex manual.tex
	pdflatex manual.tex

apply_Tstar.pdf: apply_Tstar.tex
	pdflatex apply_Tstar.tex
	pdflatex apply_Tstar.tex

%.pdf: %.eps
	epstopdf $< --outfile=$@

clean:
	rm -f *.aux *.log *.tex~