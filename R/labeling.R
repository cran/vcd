#################################################################
## labeling

pexpand <- function(par, len, default_value, default_names, choices = NULL) {
  if (!is.null(choices))
    par <- as.character(sapply(par, match.arg, choices))
  nam <- names(par)
  if (is.null(nam))
    default_value <- par
  else if (length(nam[nam == ""])) {
    default_value <- par[nam == ""]
    nam <- nam[nam != ""]
  }
  ret <- rep(default_value, length.out = len)
  if (!is.null(nam)) {
    names(ret) <- default_names
    ret[nam] <- par[nam]
    ret
  } 
  ret
}

labeling_list <- function(gp = gpar(),
                          just = "left",
                          pos = "left",
                          lsep = ": ", sep = " ",
                          offset = unit(c(2, 2), "lines"),
                          varnames = TRUE,
                          cols = 2,
                          ...) 
  function(d, split_vertical, condvars) {
    if (is.table(d))
      d <- dimnames(d)
    ld <- length(d)
    labeling_text(labels = FALSE, varnames = varnames)(d, split_vertical, condvars)
    seekViewport("marginBottom")
    pos <- unit(switch(pos, left = 0, center = 0.5, 1) / cols, "npc")
    ind <- split(seq(ld), rep.int(seq(cols), ceiling(ld / cols))[seq(ld)])
    
    for (i in seq(along = ind))
      grid.text(x = offset[1] + pos + unit((i - 1) / cols, "npc"),
                y = unit(1, "npc") - offset[2],
                paste(names(d[ind[[i]]]),
                      sapply(d[ind[[i]]], paste, collapse = sep),
                      sep = lsep,
                      collapse = "\n"
                      ),
                just = c(just, "top"),
                gp = gp
                )
  }
class(labeling_list) <- "panel_generator"

labeling_conditional <- function(...)
  function (d, split_vertical, condvars) {
    if (is.table(d))
      d <- dimnames(d)
    v <- rep.int(TRUE, length(d))
    v[seq(condvars)] <- FALSE
    labeling_text(labels = !v, ...)(d, split_vertical, condvars)
    labeling_cells(labels = v, ...)(d, split_vertical, condvars)
  }
class(labeling_conditional) <- "panel_generator"

labeling_cells <- function(labels = TRUE, varnames = TRUE,
                         abbreviate_labels = FALSE, abbreviate_varnames = FALSE,
                         gp = gpar(), lsep = ": ", lcollapse = "\n",
                         just = "center", pos = "center", rot = 0,
                         margin = unit(0.5, "lines"), clip_cells = TRUE,
                         text = NULL, ...)
  function(d, split_vertical, condvars) {
    if (is.table(d))
      d <- dimnames(d)
    dn <- names(d)
    ld <- length(d)

    ## expand parameters
    if (length(pos) < 2) pos <- c(pos, pos)
    labels <- pexpand(labels, ld, TRUE, dn)
    varnames <- pexpand(varnames, ld, TRUE, dn)
    abbreviate_labels <- pexpand(abbreviate_labels, ld, FALSE, dn)
    abbreviate_varnames <- pexpand(abbreviate_varnames, ld, FALSE, dn)

    ## margin
    if (!is.unit(margin))
      margin <- unit(margin, "lines")
    
    prvars <- ifelse(abbreviate_varnames,
                     sapply(seq(along = dn),
                            function(i) abbreviate(dn[i], abbreviate_varnames[i])),
                     dn)
    prvars <- ifelse(varnames, paste(prvars, lsep, sep = ""), "")

    ## draw labels
    split <- function(vind = 1, labs = c()) {
      n <- d[[vind]]
      for (labind in seq(along = n)) {
        lab <- c(labs, n[labind])
        names(lab) <- names(d)[1:vind]
        mlab <- paste("cell", paste(dn[1:vind], lab, sep = ".", collapse = ".."),
                      sep = "..")

        if (vind < ld)
          split(vind + 1, lab)
        else {
          seekViewport(mlab)
          pushViewport(viewport(width = max(unit(0, "npc"), unit(1, "npc") - 2 * margin),
                                height = unit(1, "npc") - 2 * margin,
                                clip = clip_cells))
          txt <- if (!is.null(text)) {
            lab <- lab[names(dimnames(text))]
            do.call("[", c(list(text), as.list(lab)))
          } else {
            prlab <- ifelse(abbreviate_labels,
                            sapply(seq(along = lab),
                                   function(i) abbreviate(lab[i], abbreviate_labels[i])),
                            lab)
            prlab <- prlab[labels[1:ld]]
            paste(prvars[labels[1:ld]], prlab, sep = "", collapse = lcollapse) 
          } 
          
          grid.text(if(!is.na(txt)) txt,
                    x = switch(pos[1], left =, top = 0, center = 0.5, 1),
                    y = switch(pos[2], left =, top = 1, center = 0.5, 0),
                    gp = gp, just = just, rot = rot)
          popViewport()
        }
      }
    }
    split()
    
}
class(labeling_cells) <- "panel_generator"

labeling_text <- function(labels = TRUE, varnames = labels,
                          tl_labels = NULL, tl_varnames = NULL, 
                          gp_labels = gpar(fontsize = 12),
                          gp_varnames = gpar(fontsize = 12, fontface = 2),
                          rot_labels = c(0, 90, 0, 90),
                          rot_varnames = c(0, 90, 0, 90),
                          pos_labels = "center", pos_varnames = "center",
                          just_labels = "center", just_varnames = pos_varnames,
                          boxes = FALSE, fill_boxes = NULL,
                          offset = c(0, 0, 0, 0),
                          
                          labbl_varnames = NULL,
                          labels_varnames = FALSE, sep = ": ",
                          
                          abbreviate = FALSE, rep = TRUE,
                          clip = FALSE, ...
                          )
  function(d, split_vertical, condvars) {
    if (is.table(d))
      d <- dimnames(d)
    dn <- names(d)
    ld <- length(d)

    ## expand parameters
    clip <- pexpand(clip, ld, TRUE, dn)
    labels <- pexpand(labels, ld, TRUE, dn)
    labels_varnames <- pexpand(labels_varnames, ld, FALSE, dn)
    pos_labels <- pexpand(pos_labels, 4, "center", c("top", "right", "bottom", "left"), c("left", "center", "right"))
    just_labels <- pexpand(just_labels, 4, "center", c("top", "right", "bottom", "left"), c("left", "center", "right"))
    offset <- if (!is.unit(offset))
      unit(pexpand(offset, 4, rep.int(0, 4), c("top","right","bottom","left")), "lines")
    else
      unit.rep(offset, length.out = 4)

    ## tl_labels
    def <- logical()
    def[split_vertical] <- rep(c(TRUE, FALSE), length.out = sum(split_vertical))
    def[!split_vertical] <- rep(c(TRUE, FALSE), length.out = sum(!split_vertical))
    tl_labels <- if (is.null(tl_labels)) 
      def
    else
      pexpand(tl_labels, ld, def, dn)

    ## rep labels
    rep <- pexpand(rep, ld, TRUE, dn)
    printed <- lapply(d, function(i) rep.int(FALSE, length(i)))
    
    ## abbreviate
    abbreviate <- pexpand(abbreviate, ld, FALSE, dn)
    labs <- d
    for (i in seq(along = d))
      if (abbreviate[i])
        labs[[i]] <- abbreviate(labs[[i]], abbreviate[i])

    ## gp_labels
    if (inherits(gp_labels, "gpar"))
      gp_labels <- list(gp_labels)
    gp_labels <- pexpand(gp_labels, ld, gpar(fontsize = 12), dn)

    ## rot_labels: top/right/bottom/left
    rot_labels <- pexpand(rot_labels, 4, c(0, 90, 0, 90),
                          c("top", "right", "bottom", "left"))

    ## varnames
    varnames <- pexpand(varnames, ld, labels, dn)

    ## gp_varnames: top/right/bottom/left!
    if (inherits(gp_varnames, "gpar"))
      gp_varnames <- list(gp_varnames)
    gp_varnames <- pexpand(gp_varnames, 4, gpar(fontsize = 12, fontface = 2),
                           c("top", "right", "bottom", "left"))

    ## rot_varnames: top/right/bottom/left!
    rot_varnames <- pexpand(rot_varnames, 4, c(0, 90, 0, 90),
                          c("top", "right", "bottom", "left"))

    ## pos_varnames: top/right/bottom/left!
    pos_varnames <- pexpand(pos_varnames, 4, "center",
                           c("top", "right", "bottom", "left"),
                            c("left", "center", "right"))

    ## just_varnames: top/right/bottom/left!
    just_varnames <- pexpand(just_varnames, 4, pos_varnames,
                             c("top", "right", "bottom", "left"),
                             c("left", "center", "right"))

    ## tl_varnames
    if (is.null(tl_varnames) && is.null(labbl_varnames))
      tl_varnames <- tl_labels
    tl_varnames <- pexpand(tl_varnames, ld, tl_labels, dn)

    ## labbl_varnames
    if (!is.null(labbl_varnames))
      labbl_varnames <- pexpand(labbl_varnames, ld, TRUE, dn)

    ## boxes
    boxes <- pexpand(boxes, ld, FALSE, dn)

    ## fill_boxes
    if (is.null(fill_boxes))
      fill_boxes <- lapply(sapply(d, length),
                           function(i) gray(0.3 + 0.4 * rev(seq(i)) / i))
    else
      fill_boxes <- pexpand(fill_boxes, ld, "grey", dn)

    ## precompute spaces
    lsp <- tsp <- bsp <- rsp <- 0
    labsp <- rep.int(0, ld)
    for (i in seq(along = dn)[tl_labels & labels])
      labsp[i] <- if (split_vertical[i])
        tsp <- tsp + 1
      else
        lsp <- lsp - 1
    for (i in rev(seq(along = dn)[!tl_labels & labels]))
      labsp[i] <- if (split_vertical[i])
        bsp <- bsp - 1
      else
        rsp <- rsp + 1
    
    if(is.null(labbl_varnames)) {
    ## varnames in the outer margin  
      ## compute axis names
      tt <- bt <- lt <- rt <- ""
      for (i in seq(along = dn))
        if (varnames[i]) {
          if (split_vertical[i]) {
            if (tl_varnames[i])
              tt <- paste(tt, dn[i], sep = if (tt == "") "" else " / ")
            else
              bt <- paste(bt, dn[i], sep = if (bt == "") "" else " / ")
          } else {
            if (tl_varnames[i])
              lt <- paste(lt, dn[i], sep = if (lt == "") "" else " / ")
            else
              rt <- paste(rt, dn[i], sep = if (rt == "") "" else " / ")
          }
        }

      ## draw axis names
      if (tt != "")
        grid.text(tt, y = unit(1, "npc") + unit(tsp + 1, "lines") + offset[1],
                  x = switch(pos_varnames[1], left =, bottom = 0, center =, centre = 0.5, 1),
                  rot = rot_varnames[1], just = just_varnames[1], gp = gp_varnames[[1]])
      if (bt != "")
        grid.text(bt, y = unit(bsp - 1, "lines") + -1 * offset[3],
                  x = switch(pos_varnames[3], left =, bottom = 0, center =, centre = 0.5, 1),
                  rot = rot_varnames[3], just = just_varnames[3], gp = gp_varnames[[3]])
      if (lt != "")
        grid.text(lt, x = unit(lsp - 1, "lines") + -1 * offset[4],
                  y = switch(pos_varnames[4], left =, bottom = 0, center =, centre = 0.5, 1),
                  rot = rot_varnames[4], just = just_varnames[4], gp = gp_varnames[[4]])
      if (rt != "")
        grid.text(rt, x = unit(1, "npc") + unit(rsp + 1, "lines") + offset[2],
                  y = switch(pos_varnames[2], left =, bottom = 0, center =, centre = 0.5, 1),
                  rot = rot_varnames[2], just = just_varnames[2], gp = gp_varnames[[2]])
    } else {
    ## varnames beneath labels
      for (i in seq(along = dn))
        if (varnames[i]) {
          if (split_vertical[i]) {
            if (tl_labels[i]) {
              if (labbl_varnames[i]) {
                grid.text(dn[i],
                          y = unit(1, "npc") + unit(1 + tsp - labsp[i], "lines") + offset[1],
                          x = unit(-0.5, "lines"),
                          just = "right", gp = gpar(fontface = 2))
              } else {
                grid.text(dn[i], y = unit(1, "npc") + unit(1 + tsp - labsp[i], "lines") + offset[1],
                          x = unit(1, "npc") + unit(0.5, "lines"),
                          just = "left", gp = gpar(fontface = 2))
              }
            } else {
              if (labbl_varnames[i]) {
                grid.text(dn[i], y = unit(labsp[i], "lines") + -1 * offset[3],
                          x = unit(-0.5, "lines"), just = "right",
                          gp = gpar(fontface = 2))
              } else {
                grid.text(dn[i], y = unit(labsp[i], "lines") + -1 * offset[3],
                          x = unit(1, "npc") + unit(0.5, "lines"),
                          just = "left", gp = gpar(fontface = 2))
              }
            }
          } else {
            if (tl_labels[i]) {
              if (labbl_varnames[i]) {
                grid.text(dn[i], x = unit(lsp - 1 - labsp[i], "lines") + -1 * offset[4],
                          y = unit(-0.5, "lines"), just = "right", rot = 90,
                          gp = gpar(fontface = 2))
              } else {
                grid.text(dn[i], x = unit(lsp - 1 - labsp[i], "lines") + -1 * offset[4],
                          y = unit(1, "npc") + unit(0.5, "lines"),
                          just = "left", rot = 90, gp = gpar(fontface = 2))
              }
            } else {
              if (labbl_varnames[i]) {
                grid.text(dn[i], x = unit(1, "npc") + unit(labsp[i], "lines") + offset[2],
                          y = unit(-0.5, "lines"),
                          just = "right", rot = 90, gp = gpar(fontface = 2))
              } else {
                grid.text(dn[i], x = unit(1, "npc") + unit(labsp[i], "lines") + offset[2],
                          y = unit(1, "npc") + unit(0.5, "lines"),
                          just = "left", rot = 90, gp = gpar(fontface = 2))
              }
            }
          }
        }
    }

    ## draw labels
    split <- function(vind = 1, root = "cell",
                      left = TRUE, right = TRUE, top = TRUE, bottom = TRUE) {
      n <- d[[vind]]
      vl <- length(n)
      sp <- split_vertical[vind]
      labseq <- seq(along = n)
      if (!sp) labseq <- rev(labseq)
      
      for (labind in labseq) {
        mlab <- paste(root, "", dn[vind], n[labind], sep = ".")
        if (labels[vind] && (rep[vind] || !printed[[vind]][labind])) {
          lab <- labs[[vind]][labind]
          if (labels_varnames[vind])
            lab <- paste(dn[vind], lab, sep = sep)
          if (sp) {
            if (tl_labels[vind]) {
              if (top) {
                seekViewport(mlab)
                if (clip[vind])
                  pushViewport(viewport(height = unit(1, "npc") + 2 * offset[1] +
                                        unit(2 * (2 + tsp - labsp[vind]), "lines"),
                                        clip = "on"))
                if (boxes[vind])
                  grid.rect(height = unit(0.8, "lines"),
                            y = unit(1, "npc") + offset[1] +
                                unit(1 + tsp - labsp[vind] - (2 + as.numeric(offset[1]) + tsp - labsp[vind]) * clip[vind], "lines"),
                            gp = gpar(fill = fill_boxes[[vind]][labind]))
                grid.text(lab,
                          y = unit(1, "npc") + offset[1] +
                              unit(1 + tsp - labsp[vind] - (2 + as.numeric(offset[1]) + tsp - labsp[vind]) * clip[vind], "lines"),
                          x = unit(0.15 * switch(pos_labels[1], left =, bottom = 1, center =, centre = 0, -1) * boxes[vind], "lines") +
                              unit(switch(pos_labels[1], left =, bottom = 0, center =, centre = 0.5, 1), "npc"),
                          rot = rot_labels[1], just = just_labels[1],
                          gp = gp_labels[[vind]])
                if (clip[vind]) popViewport()
                printed[[vind]][labind] <<- TRUE
              }
            } else {
              if (bottom) {
                seekViewport(mlab)
                if (clip[vind])
                  pushViewport(viewport(height = unit(1, "npc") + 2 * offset[3] + 
                                        unit(2 * (1 + abs(labsp[vind])), "lines"),
                                        clip = "on"))
###
                if (boxes[vind])
                  grid.rect(height = unit(0.8, "lines"),
                            y = -1 * offset[3] + unit(labsp[vind] + (1 + as.numeric(offset[3]) + abs(labsp[vind])) * clip[vind], "lines"),
                            gp = gpar(fill = fill_boxes[[vind]][labind]))

                grid.text(lab,
                          y = -1 * offset[3] + unit(labsp[vind] + (1 + as.numeric(offset[3]) + abs(labsp[vind])) * clip[vind], "lines"),
                          x = unit(0.15 * switch(pos_labels[3], left =, bottom = 1, center =, centre = 0, -1) * boxes[vind], "lines") +
                              unit(switch(pos_labels[3], left =, bottom = 0, center =, centre = 0.5, 1), "npc"),
                          rot = rot_labels[3], just = just_labels[3],
                          gp = gp_labels[[vind]])
                if (clip[vind]) popViewport()
                printed[[vind]][labind] <<- TRUE
              }
            }
          } else {
            if (tl_labels[vind]) {
              if (left) {
                seekViewport(mlab)
                if (clip[vind])
                  pushViewport(viewport(width = unit(1, "npc") + 2 * offset[4] + 
                                        unit(2 * (2 - lsp + labsp[vind]), "lines"),
                                        clip = "on"))
                if (boxes[vind])
                  grid.rect(width = unit(0.8, "lines"),
                            x = -1 * offset[4] + unit(lsp - 1 - labsp[vind] + (2 - lsp + as.numeric(offset[4]) + labsp[vind]) * clip[vind], "lines"),
                            gp = gpar(fill = fill_boxes[[vind]][labind]))

                grid.text(lab,
                          x = -1 * offset[4] + unit(lsp - 1 - labsp[vind] + (2 - lsp + as.numeric(offset[4]) + labsp[vind]) * clip[vind], "lines"),
                          y = unit(0.15 * switch(pos_labels[4], left =, bottom = 1, centre = 0, -1) * boxes[vind], "lines") +
                          unit(switch(pos_labels[4], left =, bottom = 0, center =, centre = 0.5, 1), "npc"),
                          rot = rot_labels[4], just = just_labels[4],
                          gp = gp_labels[[vind]])

                if (clip[vind]) popViewport()
                printed[[vind]][labind] <<- TRUE
              }
            } else {
              if (right) {
                seekViewport(mlab)
                if (clip[vind])
                  pushViewport(viewport(width = unit(1, "npc") + 2 * offset[2] +
                                        unit(2 * (1 + abs(labsp[vind])), "lines"),
                                        clip = "on"))
                if (boxes[vind])
                  grid.rect(width = unit(0.8, "lines"),
                            x = offset[2] + unit(1, "npc") +
                                unit(labsp[vind] - (1 + as.numeric(offset[2]) + abs(labsp[vind])) * clip[vind], "lines"),
                            gp = gpar(fill = fill_boxes[[vind]][labind]))
                grid.text(lab,
                          x = offset[2] + unit(1, "npc") + unit(0.1, "lines") +
                              unit(labsp[vind] - (1 + as.numeric(offset[2]) + abs(labsp[vind])) * clip[vind], "lines"),
                          y = unit(0.15 * switch(pos_labels[2], left =, bottom = 1, center =, centre = 0, -1) * boxes[vind], "lines") +
                              unit(switch(pos_labels[2], left =, bottom = 0, center =, centre = 0.5, 1), "npc"),
                          rot = rot_labels[2], just = just_labels[2],
                          gp = gp_labels[[vind]])

                if (clip[vind]) popViewport()
                printed[[vind]][labind] <<- TRUE
              }
            }
          }
        }
        
        if (vind < ld) Recall(vind + 1, mlab,
                              if (sp) left && labind == 1 else left,
                              if (sp) right && labind == vl else right,
                              if (!sp) top && labind == 1 else top,
                              if (!sp) bottom && labind == vl else bottom)
      }
    }
    split()
    
  }
class(labeling_text) <- "panel_generator"

labeling_doubledecker <- function(lab_pos = c("bottom", "top"), ...) {
  lab_pos <- match.arg(lab_pos)
  function(d, split_vertical, condvars) {
    if (is.table(d))
      d <- dimnames(d)
    labeling_text(boxes = c(rep.int(TRUE, length(d) - 1), FALSE),
                  clip = c(rep.int(TRUE, length(d) - 1), FALSE),
                  labbl_varnames = FALSE,
                  rot_labels = rep.int(0, 4),
                  pos_labels = c("left", "center", "left", "center"),
                  just_labels = c("left", "left", "left", "center"),
                  varnames = c(c(rep.int(TRUE, length(d) - 1), FALSE)),
                  offset = c(0, -0.6, 0, 0),
                  tl_labels = c(rep.int(lab_pos== "top", length(d) - 1), FALSE)
                  )(d, split_vertical, condvars)
    seekViewport("marginRight")
    grid.text(names(d)[length(d)],
              x = unit(0.5, "lines"), y = unit(1, "npc"), just = c("left","top"),
              gp = gpar(fontface = 2))
  }
}
class(labeling_doubledecker) <- "panel_generator"

labeling_left <- function(tl_labels = TRUE, clip = TRUE, pos_varnames = "left",
                        pos_labels = "left", just_labels = "left", ...)
  labeling_text(tl_labels = tl_labels, clip = clip, pos_varnames = pos_varnames,
              pos_labels = pos_labels, just_labels = just_labels, ...)
class(labeling_left) <- "panel_generator"

labeling_cboxed <- function(tl_labels = TRUE, boxes = TRUE, clip = TRUE, pos_labels = "center", ...)
  labeling_text(tl_labels = tl_labels, boxes = boxes, clip = clip, pos_labels = pos_labels, ...)
class(labeling_cboxed) <- "panel_generator"

labeling_lboxed <- function(tl_labels = FALSE, boxes = TRUE, clip = TRUE, pos_labels = "left", just_labels = "left", labbl_varnames = FALSE, ...)
  labeling_text(tl_labels = tl_labels, boxes = boxes, clip = clip, pos_labels = pos_labels, labbl_varnames = labbl_varnames, just_labels = just_labels, ...)
class(labeling_lboxed) <- "panel_generator"
