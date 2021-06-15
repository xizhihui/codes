#' ggaxispoint
#' 
#' Add plots(like point) to the axis text. 
#' This code is modified from \code{\link{ggtext::element_markdown} and \code{\link{gridtext::richtext_grob}}.
#' 
NULL

#' Theme element that enables addional plots(points, polygon, etc.).
#' 
#' Theme element that enables addional plots(points, polygon, etc.).
#' 
#' @param family Font family
#' @param face Font face
#' @param size Font size
#' @param colour,color Text color
#' @param fill Fill color of the enclosing box
#' @param linetype Line type of the enclosing box (like `lty` in base R)
#' @param linewidth Line width of the enclosing box (measured in mm, just like `size` in
#'   [ggplot2::element_line()]).
#' @param hjust Horizontal justification
#' @param vjust Vertical justification
#' @param halign Horizontal justification
#' @param valign Vertical justification
#' @param lineheight Line height
#' @param padding,margin Padding and margins around the text box.
#'   See [gridtext::richtext_grob()] for details.
#' @param r Unit value specifying the corner radius of the box 
#' @param angle Angle (in degrees)
#' @param align_widths,align_heights Should multiple elements be aligned by their
#'   widths or height? See [gridtext::richtext_grob()] for details.
#' @param rotate_margins Should margins get rotated in frame with rotated text?
#'   If `TRUE`, the margins are applied relative to the text direction. If `FALSE`,
#'   the margins are applied relative to the plot direction, i.e., the top margin,
#'   for example, is always placed above the text label, regardless of the direction
#'   in which the text runs. The default is `FALSE`, which mimics the behavior of 
#'   `element_text()`.
#' @param debug Draw a debugging box around each label
#' @param inherit.blank See [ggplot2::margin()] for details.
#' @return A ggplot2 theme element that can be used inside a [ggplot2::theme()]
#'   call.
#' @seealso [gridtext::richtext_grob()], [element_textbox()], [geom_richtext()]
#' @export
element_axispoint <- function(
    family = NULL, face = NULL, size = NULL, colour = NULL, fill = NULL,
    linetype = NULL, linewidth = NULL,
    hjust = NULL, vjust = NULL, halign = NULL, valign = NULL, angle = NULL,
    lineheight = NULL, margin = NULL, padding = NULL, r = NULL, 
    color = NULL, box.color = NULL, align_widths = NULL,
    align_heights = NULL, rotate_margins = NULL,
    debug = FALSE, inherit.blank = FALSE,
    ap.position = c("inner", "outer"), ap.size = 1, ap.color = "green"
) {
    if (!is.null(color)) colour <- color
    ap.position <- match.arg(ap.position)
    
    structure(
        list(
            family = family, face = face, size = size, colour = colour, fill = fill,
            linetype = linetype, linewidth = linewidth,
            hjust = hjust, vjust = vjust, halign = halign, valign = valign, angle = angle,
            lineheight = lineheight, margin = margin, padding = padding, r = r,
            color = color, align_widths = align_widths,
            align_heights = align_heights, rotate_margins = rotate_margins,
            debug = FALSE, inherit.blank = inherit.blank,
            ap.position = ap.position, ap.size = ap.size, ap.color = ap.color
        ),
        class = c("element_axispoint", "element_text", "element")
    )
}

element_grob.element_axispoint <- function(
    element, label = "", x = NULL, y = NULL,
    family = NULL, face = NULL, colour = NULL, size = NULL,
    hjust = NULL, vjust = NULL, angle = NULL, lineheight = NULL,
    margin = NULL, margin_x = FALSE, margin_y = FALSE,
    ap.position = NULL, ap.size = NULL, ap.color = NULL,
    ...
) {
    if (is.null(label)) return(ggplot2::zeroGrob())
    
    `%||%` <- ggplot2:::`%||%`
    hj <- hjust %||% element$hjust
    vj <- hjust %||% element$vjust
    halign <- element$halign %||% hj
    valign <- element$valign %||% vj
    padding <- element$padding %||% ggplot2::margin(0, 0, 0, 0)
    margin <- margin %||% element$margin %||% ggplot2::margin(0, 0, 0, 0)
    angle <- angle %||% element$angle %||% 0
    r <- element$r %||% unit(0, "pt")
    ap.position <- ap.position %||% element$ap.position
    ap.size <- ap.size %||% element$ap.size
    ap.color <- ap.color %||% element$ap.color
    align_widths <- isTRUE(element$align_widths)
    align_heights <- isTRUE(element$align_heights)
    
    # We rotate the justifiation values to obtain the correct x and y reference point,
    # since box_hjust and box_vjust are applied relative to the rotated text frame in richtext_grob
    # for me: Is this required for axispoint? To be determined.
    just <- ggtext:::rotate_just(angle, hj, vj)
    
    n <- max(length(x), length(y), 1)
    x <- x %||% unit(rep(just$hjust, n), "npc")
    y <- y %||% unit(rep(just$vjust, n), "npc")
    
    # The gp settings can override element_gp.
    # for me: how to use element_gp?
    txtgp <- grid::gpar(
        fontsize = size %||% element$size,
        col = colour %||% element$colour,
        fontfamily = family %||% element$family,
        fontface = face %||% element$face,
        lineheight = lineheight %||% element$lineheight
    )
    
    grob <- make_axispoint_grob(
        label, x = x, y = y, hjust = hj, vjust = vj,
        halign = halign, valign = valign, rot = angle,
        padding = padding, margin = unit(c(0, 0, 0, 0), "pt"), r = r,
        align_widths = align_widths, align_heights = align_heights,
        gp = txtgp, debug = element$debug,
        ap.position = ap.position,
        ap.size = ap.size,
        ap.color = ap.color
    )
    grob <- ggtext:::add_margins(grob, margin, margin_x, margin_y, debug = element$debug)
    grob
}

make_axispoint_grob <- function(
    text, x = unit(0.5, "npc"), y = unit(0.5, "npc"),
    hjust = 0.5, vjust = 0.5, halign = hjust, valign = vjust,
    rot = 0, default.units = "npc",
    margin = unit(c(0, 0, 0, 0), "pt"), padding = unit(c(0, 0, 0, 0), "pt"),
    r = unit(0, "pt"), align_widths = FALSE, align_heights = FALSE,
    name = NULL, gp = grid::gpar(), vp = NULL, debug = FALSE,
    ap.position = NULL,
    ap.size = NULL,
    ap.color = NULL
) {
    #browser()
    if (!grid::is.unit(x)) x <- unit(x, default.units)
    if (!grid::is.unit(y)) y <- unit(y, default.units)
    text <- as.character(text)
    text <- ifelse(is.na(text), "", text)
    
    # make sure text, x, and y have the same length
    n <- unique(length(text), length(x), length(y))
    if (length(n) > 1) {
        stop("Arguments `text`, `x`, and `y` must have the same length.", call. = FALSE)
    }
    if (length(ap.color) == 1) {
        ap.color <- rep(ap.color, length(text))
    }
    if (length(ap.size) == 1) {
        ap.size <- rep(ap.size, length(text))
    }
    
    if (length(unique(x)) == 1) {
        current_axis <- "y"
    } else {
        current_axis <- "x"
    }
    
    gp_list <- gridtext:::recycle_gpar(gp, n)
    
    text_x_list <- gridtext:::unit_to_list(x)
    text_y_list <- gridtext:::unit_to_list(y)
    
    fix_position <- function(raw, offset, position) {
        if (position == "outer") {
            position <- list(text = raw, point = raw - offset * 1.5)
        } else {
            position <- list(text = raw - offset * 1.5, point = raw - offset)
        }
        position
    }
    if (current_axis == "x") {
        #offset <- lapply(text_grob, grid::heightDetails)
        offset <- lapply(text, grid::stringHeight)
        position_x_text <- text_x_list
        position_x_point <- text_x_list
        position_y <- mapply(
            fix_position,
            raw = text_y_list,
            offset = offset,
            position = ap.position,
            SIMPLIFY = FALSE
        )
        position_y_text <- lapply(position_y, function(x) x$text)
        position_y_point <- lapply(position_y, function(x) x$point)
    } else {
        #offset <- lapply(text_grob, grid::widthDetails)
        offset <- lapply(text, grid::stringWidth)
        position_y_text <- text_y_list
        position_y_point <- text_y_list
        position_x <- mapply(
            fix_position,
            raw = text_x_list,
            offset = offset,
            position = ap.position,
            SIMPLIFY = FALSE
        )
        position_x_text <- lapply(position_x, function(x) x$text)
        position_x_point <- lapply(position_x, function(x) x$point)
    }
    
    text_grob <- mapply(
        grid::textGrob,
        label = text,
        x = position_x_text,
        y = position_y_text,
        hjust = hjust,
        vjust = vjust,
        rot = rot,
        default.units = default.units,
        gp = gp_list,
        SIMPLIFY = FALSE
    )
    text_grob_width <- lapply(text_grob, grid::grobWidth)
    text_grob_height <- lapply(text_grob, grid::grobHeight)
    
    points_grob <- mapply(
        grid::pointsGrob,
        x = position_x_point,
        y = position_y_point,
        pch = 19,
        gp = lapply(ap.color, function(col) grid::gpar(col = col)),
        size = lapply(ap.size, grid::unit, units = "char"),
        default.units = default.units,
        SIMPLIFY = FALSE
    )
    points_grob_width <- lapply(points_grob, grid::grobWidth)
    points_grob_height <- lapply(points_grob, grid::grobHeight)
    
    # below is not worked
    if (ap.position == "inner") {
        children <- c(text_grob, points_grob)
    } else {
        children <- c(points_grob, text_grob)
    }
    children <- do.call(grid::gList, children)
    fix_width_height <- function(x, y, axis, wh) {
        xwidth <- axis == "x" && wh == "width"
        #xheight <- axis == "x" && wh == "height"
        #ywidth <- axis == "y" && wh == "width"
        yheight <- axis == "y" && wh == "height"
        
        if (xwidth || yheight) {
            max(x, y)
        } else {
            x + 2 * y
        }
    }
    children_width <- mapply(
        fix_width_height,
        x = text_grob_width,
        y = points_grob_width,
        MoreArgs = list(
            axis = current_axis,
            wh = "width"
        ),
        SIMPLIFY = FALSE
    )
    children_width <- do.call(max, children_width)
    children_height <- mapply(
        fix_width_height,
        x = text_grob_height,
        y = points_grob_height,
        MoreArgs = list(
            axis = current_axis,
            wh = "height"    
        ),
        SIMPLIFY = FALSE
    )
    children_height <- do.call(max, children_height)
    
    grid::gTree(
        gp = gp,
        vp = vp,
        name = name,
        wext = children_width,
        hext = children_height,
        debug = debug,
        children = children,
        cl = "axispoint_grob"
    )
}

heightDetails.axispoint_grob <- function(x) {
    do.call(max, lapply(x$children, grid::grobHeight)) * 2
    #x$hext
}
widthDetails.axispoint_grob <- function(x) {
    do.call(max, lapply(x$children, grid::grobWidth)) * 2
    #x$wext
}
