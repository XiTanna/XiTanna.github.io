---
layout: archive
title: "Submitted Papers"
permalink: /submitted/
author_profile: true
---

{% if author.googlescholar %}
  
{% endif %}

{% include base_path %}

{% for post in site.submitted reversed %}
  {% include archive-single.html %}
{% endfor %}