sim.choice <- function() {
  # question text
  questions <- c(
    snps = "Do you have SNP data?",
    non.diploid = "Is your data other than diploid?",
    marker.num = "Do you want to simulate many markers?",
    pop.size = "Do you have large population sizes?",
    complex.hist = "Do you have a complex history to simulate?",
    deep.time = "Are you looking at deep time frames?",
    demography = "Do you want to include demography?",
    management = "Does your question involve management decisions?",
    completion.time = "Do you need a short completion time?",
    computer = "Do you have large computer capacity?"
  )
  # default responses
  responses <- c(F, F, F, F, F, F, F, F, F, F)
  # response weights
  forward.wts <- c(0, 0, 0.3, 0.2, 0, 0.2, 1, 1, 0.2, 0.3)
  names(responses) <- names(forward.wts) <- names(questions)

  # loop through each question
  for(this.quest in names(questions)) {
    # create prompt with default
    prompt.def <- if(responses[this.quest]) "Y/n" else "y/N"
    prompt <- paste(questions[this.quest], " (", prompt.def, ") ", sep = "")
    # keep asking question until appropriate response received
    ans <- "empty"
    while(!ans %in% c("y", "n", "")) {
      ans <- tolower(readline(prompt))
    }
    # if not default, change response
    if(ans != "") responses[this.quest] <- ans == "y"
  }
  cat("\n")

  # find reasons that forward-time models are excluded
  fwd.excl <- forward.wts == 0 & responses
  if(any(fwd.excl)) {
    reasons <- paste(names(questions)[fwd.excl], collapse = ", ")
    cat("Forward-time simulations are excluded because: ", reasons, "\n")
  }

  # find reasons that forward-time models are required
  fwd.req <- forward.wts == 1 & responses
  if(any(fwd.req)) {
    reasons <- paste(names(questions)[fwd.req], collapse = ", ")
    cat("Forward-time simulations are required because: ", reasons, "\n")
  }

  # get relative 'score' for each model
  fwd.score <- sum(forward.wts * responses) / length(responses)
  cat("Coalescent score: ", 1 - fwd.score, "\n")
  cat("Forward-time score: ", fwd.score, "\n")

  ans <- "empty"
  prompt <- "Choose a simulator: (c)oalescent or (f)orward-time: "
  while(!ans %in% c("c", "f")) {
    ans <- tolower(readline(prompt))
  }
}
