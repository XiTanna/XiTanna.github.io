---
layout: archive
permalink: /cv/
title: "Curriculum Vitae"
---

{% include base_path %}

{% for post in site.cv reversed %}
  {% include archive-single-cv.html %}
{% endfor %}
