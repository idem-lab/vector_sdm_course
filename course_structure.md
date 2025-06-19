# Course structure:

3 day course by Nick Golding and Gerry Ryan

**Day 1:**

AM:

- Welcome / icebreaker (ROpenSci activity - stand on a line thing)
- Introductory concepts and discussion. Slides and whiteboard.

PM:

- Simulate data
  - Abundance and relative abundance
  - bias
  - PA data from planned surveys (random, biased, abundance-biased)
  - understand how to simulate presence/absence from abundance
  - understand how to calculate probability of presence from average abundance
    - `simulate_prob_presence.R`
  - PO data from presence and bias process
  - `prepare_raster_data.R` - run through but encourage students to download files from figshare as travel and bioclim files are large.
  - `simulate_data.R` - students run through alongside instructors

**Day 2:**

AM:

- Modelling
  - Logistic regression on random PA data
  - Logistic regression on presence-only with random background
  - `models.R` - students run through alongside instructors
- Theory: link functions. Whiteboard and code.
  - `link_demo.R`


PM:

- Modelling
  - Maxent presence only with random background points
  - Maxent presence only with random bg and bias layer offset
  - Students can run `models.R` alongside instructors
  - Students explore other PA data or other covariates if happy
- Theory: target-group background and bias cancellation

**Day 3:**

AM:

- Modelling
  - Fithian PA-PO-bg model
  - `models.R`!
  - continue explore alternatives from existing model set
- Theory: Fithian model
- Discussion: other topics in SDMs

PM:

- Modelling
  - own data and models
  - continue explore alternatives from existing model set
- Discussion
  - Papers
  - Own models and data
  
